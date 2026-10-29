## Figure 4 defense systems
## Format Defensefinder output for good E lenta genomes for Fadima
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(patchwork)
library(ggbeeswarm)
theme_set(theme_classic2())
defensefinder_files <- list.files("defensefinder_out/", pattern = "systems.tsv")
good_genomes <- readRDS("final_qc_metadata.rds")
good_genomes <- good_genomes[grepl("lenta", species)]
defensefinder_files <- defensefinder_files[gsub("_defense_finder.*", "", defensefinder_files) %in% good_genomes$genome]
df_data <- rbindlist(lapply(defensefinder_files, function(x){
  foo <- fread(paste0("defensefinder_out/", x))
  foo[,genome:=gsub("_defense_finder.*", "", x)]
  return(foo)
}))
fwrite(df_data, file = "defense_finder_systems_compiled_elenta.tsv", sep = "\t")

df_data <- merge(df_data, good_genomes[,list(genome, Clust2, IsMAG)], by = "genome", all.x = T)

df_data_wide <- dcast(df_data, genome+Clust2~type, value.var = "genes_count", 
                      fun.aggregate = function(x){ sum(as.numeric(x))})
df_data_all <- melt(df_data_wide, id.var = c("genome", "Clust2"), variable.name = "type")
df_data_all[type %in% c("Wadjet", "RM", "MazEF", "CasFinder", "AbiO", "AbiE", "Abi2"),wilcox.test(value[Clust2 == 1], value[Clust2 == 2])$p.value, by=type]
df_data_all[,Present:=ifelse(value > 0 , 1, 0)]

df_prev <- df_data[,length(unique(genome)), by = list(subtype, type, Clust2, IsMAG)]
df_prev <- dcast(df_prev, subtype+type~Clust2+IsMAG, value.var = "V1", fill = 0) %>% 
  melt(id.var = c("subtype", "type"), variable.factor = F)
df_prev[,c("Clust2", "IsMAG"):=tstrsplit(variable, split = "_")]
tot_counts <- good_genomes[,length(unique(genome)), by = list(Clust2, IsMAG)]
setnames(tot_counts, "V1", "TotGenomes")
df_prev[,IsMAG:=as.logical(IsMAG)]
df_prev <- merge(df_prev, tot_counts)
df_prev[,FracGenomes:=value/TotGenomes]
df_prev[,Group:=paste0(Clust2, "_", ifelse(IsMAG, "MAGs", "Isolates"))]

type_order <- df_data[,length(unique(genome)), by = type][order(V1, decreasing = T), type]

df_prev3 <- df_data[,length(unique(genome)), by = list(type, Clust2)] %>% 
  dcast(type~Clust2, value.var = "V1", fill = 0)
df_prev3[,Prev1:=`1`/149]
df_prev3[,Prev2:=`2`/135]
df_prev3[,Log2FC:=log2(Prev2/Prev1+ 1e-2)]
df_prev3[,OverallPrev:=(`1`+`2`)/284]
df_prev3[,type:=factor(type, levels = rev(type_order))]
variable_types <- df_data_all[,sum(Present), by = type][V1 > 3 & V1 < 275, type]
type_stats <- df_data_all %>% 
  merge(good_genomes[,list(genome, IsMAG, Genome_Size)]) %>% 
  filter(type %in% variable_types) %>% 
  group_by(type, IsMAG) %>% 
  summarize(broom::tidy(fisher.test(Present, Clust2))) %>% 
  mutate(p.adj = p.adjust(p.value, method = "BH")) 
type_stats %>% filter(p.adj < 0.1)
type_stats <- df_data_all %>% 
  merge(good_genomes[,list(genome, IsMAG, Genome_Size)]) %>% 
  filter(type %in% variable_types) %>% 
  group_by(type) %>% 
  summarize(broom::tidy(glm(Present~Clust2+IsMAG+Genome_Size, family = "binomial"))) %>% 
  group_by(term) %>% 
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  ungroup()
type_stats %>% filter(p.adj < 0.1 & !grepl("Intercept", term))
signif_types <- type_stats %>% filter(p.adj < 0.1 & grepl("Clust2", term)) %>% pull(type)
type_glm_est <- type_stats %>% filter(term == "Clust22")

df_prev3[,SignifColor:=case_when(
  type %in% signif_types & Log2FC > 0 ~ "Benriched", 
  type %in% signif_types & Log2FC < 0 ~ "B_depleted", 
  TRUE ~ "Not signif"
)]
df_prev3 <- merge(df_prev3, select(type_glm_est, type, estimate, p.adj), by = "type", all.x = T)

logfc_plot <- ggplot(df_prev3[OverallPrev > 0.01], aes(x = type, y = Log2FC, fill = SignifColor)) + 
  geom_bar(stat = "identity", aes(color = factor(ifelse(SignifColor == "Not signif", 0, 1)))) + coord_flip() + 
  scale_fill_manual(values = c(brewer.pal(3, "Paired")[1:2], "gray80"), guide = "none") + 
  ylab("Log2 enrichment/depletion in Lineage B genomes") + xlab("") +
  scale_color_manual(values = c("gray70", "black"), guide = "none")

prev_plot <- ggplot(df_prev3[OverallPrev > 0.01], aes(x = type, y = OverallPrev)) + 
  geom_bar(stat = "identity") + coord_flip() + ylab("Overall prevalence")

prev_types_plot <- prev_plot +scale_y_reverse() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  logfc_plot + theme(axis.text.y = element_text(hjust = 0), axis.ticks.y = element_blank())+ plot_layout(widths = c(1,2))

ggsave(prev_types_plot, file = "defense_systems_prev.pdf", width = 6, height = 4.8)


num_df_systems <- df_data[,length(subtype), by = list(genome, Clust2, IsMAG)]
num_df_systems <- merge(num_df_systems, good_genomes[,list(genome, Clust2, IsMAG, Genome_Size, n_genes_called, Completeness)], all.y = T)
num_df_systems[is.na(V1), V1:=0]

num_df_systems[,mean(V1), by = list(Clust2, IsMAG)]
# Diff in means between isolates
num_df_systems[IsMAG==F, mean(V1[Clust2==2])-mean(V1[Clust2==1])]
num_df_systems[IsMAG==T, mean(V1[Clust2==2])-mean(V1[Clust2==1])]
num_df_systems[,mean(V1), by = Clust2]
num_df_systems[,var(V1), by = Clust2] #Poisson regression appears appropriate (not overdispersed)

glm_res <- num_df_systems %>% 
  reframe(broom::tidy(glm(V1~factor(IsMAG)+factor(Clust2)+offset(log(Genome_Size)), family = "poisson"))) %>% 
  group_by(term) %>% 
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  ungroup()
glm_res %>% filter(p.adj < 0.1)

total_syst_plot <- ggplot(num_df_systems, aes(x = factor(Clust2), y = jitter(V1))) + 
  geom_beeswarm(method = "compactswarm", cex = 1.5,  aes(fill = factor(Clust2), 
                                                         shape = ifelse(IsMAG, "MAG", "Isolate")), alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, fill = NA) + ylab("Number of defense systems") +
  scale_x_discrete(name = "Genome group", labels = c("A", "B")) + scale_fill_brewer(palette = "Paired", guide = "none") + scale_shape_manual(values = c(21, 22), name = "") +
  annotate(label = paste0("p=", filter(glm_res, grepl("Clust2", term)) %>% pull(p.value) %>% round(digits = 4)), 
           geom = "text", y = 14, x = 1.5) + theme(legend.position = "bottom")

ggsave(total_syst_plot, width = 3.2, height = 3.8, file = "total_defense_systems_compareClust.pdf")

## Version splitting up isolates and MAGs
total_syst_plot <- ggplot(num_df_systems, 
                          aes(x = paste0(factor(Clust2), "_", ifelse(IsMAG, "MAG", "Isolate")), y = jitter(V1))) + 
  geom_beeswarm(method = "compactswarm", cex = 1.5,  aes(fill = factor(Clust2), 
                                                         shape = ifelse(IsMAG, "MAG", "Isolate")), alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, fill = NA) + ylab("Number of defense systems") +
  scale_x_discrete(name = "Genome group", labels = c("A_Isolate","A_MAG", "B_Isolate", "B_MAG")) + 
  scale_fill_brewer(palette = "Paired", guide = "none") + scale_shape_manual(values = c(21, 22), name = "") +
  annotate(label = paste0("p=", filter(glm_res, grepl("Clust2", term)) %>% pull(p.value) %>% round(digits = 4)), 
           geom = "text", y = 14, x = 2.5) + theme(legend.position = "bottom")

ggsave(total_syst_plot, width = 3.5, height = 3.8, file = "total_defense_systems_compareClust_isolates_MAGs.pdf")


total_syst_plot <- ggplot(num_df_systems, 
                          aes(x = factor(Clust2), y = jitter(V1))) + 
  geom_beeswarm(method = "compactswarm", cex = 1.5,  aes(fill = factor(Clust2), 
                                                         shape = ifelse(IsMAG, "MAG", "Isolate")), alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, fill = NA) + ylab("Number of defense systems") +
  scale_x_discrete(name = "Genome group", labels = c("A", "B")) + 
  facet_grid(~ifelse(IsMAG, "MAGs", "Isolates")) + 
  scale_fill_brewer(palette = "Paired", guide = "none") + scale_shape_manual(values = c(21, 22), name = "") +
  annotate(label = paste0("Group B p=", filter(glm_res, grepl("Clust2", term)) %>% pull(p.value) %>% round(digits = 4)), 
           geom = "text", y = 14, x = 2) + theme(legend.position = "bottom") +
  theme(strip.background = element_rect(color = "white"))
ggsave(total_syst_plot, width = 3.5, height = 3.8, file = "total_defense_systems_compareClust_isolates_MAGs.pdf")

## Test isolates alone
num_df_systems[IsMAG==F, wilcox.test(V1[Clust2==1], V1[Clust2==2])]
# p = 0.0.99

rm_counts_per_genome <- df_data[type == "RM", length(subtype), by = list(genome, Clust2, IsMAG)]
rm_counts_per_genome <- merge(rm_counts_per_genome, good_genomes[,list(genome, Clust2, IsMAG, Genome_Size, n_genes_called, Completeness)], all.y = T)
rm_counts_per_genome[is.na(V1), V1:=0]

rm_counts_per_genome[,mean(V1), by = Clust2]
rm_counts_per_genome[,var(V1), by = Clust2] #Poisson regression appears appropriate (not overdispersed)

rm_counts_per_genome[IsMAG==F, mean(V1[Clust2==2])-mean(V1[Clust2==1])]
rm_counts_per_genome[IsMAG==T, mean(V1[Clust2==2])-mean(V1[Clust2==1])]

glm_res <- rm_counts_per_genome %>% 
  reframe(broom::tidy(glm(V1~factor(IsMAG)+factor(Clust2)+offset(log(Genome_Size)), family = "poisson"))) %>% 
  group_by(term) %>% 
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  ungroup()
glm_res %>% filter(p.adj < 0.1)


total_rm_plot <- ggplot(rm_counts_per_genome, aes(x = factor(Clust2), y = jitter(V1))) + 
  geom_beeswarm(method = "compactswarm", cex = 1.5,  aes(fill = factor(Clust2), 
                                                         shape = ifelse(IsMAG, "MAG", "Isolate")), alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, fill = NA) + ylab("Number of RM systems") +
  scale_x_discrete(name = "Genome group", labels = c("A", "B")) + scale_fill_brewer(palette = "Paired", guide = "none") + scale_shape_manual(values = c(21, 22), name = "") +
  annotate(label = paste0("p=", filter(glm_res, grepl("Clust2", term)) %>% pull(p.value) %>% round(digits = 6)), 
           geom = "text", y = 5.2, x = 1.5) + theme(legend.position = "bottom")
ggsave(total_rm_plot, height = 3.8, width = 3.1, file = "total_rm_systems_clust2_boxplot.pdf")    

# Version splitting up MAGs and Isolates
total_rm_plot <- ggplot(rm_counts_per_genome, 
                          aes(x = paste0(factor(Clust2), "_", ifelse(IsMAG, "MAG", "Isolate")), y = jitter(V1))) + 
  geom_beeswarm(method = "compactswarm", cex = 1.5,  aes(fill = factor(Clust2), 
                                                         shape = ifelse(IsMAG, "MAG", "Isolate")), alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, fill = NA) + ylab("Number of RM systems") +
  scale_x_discrete(name = "Genome group", labels = c("A_Isolate","A_MAG", "B_Isolate", "B_MAG")) + 
  scale_fill_brewer(palette = "Paired", guide = "none") + scale_shape_manual(values = c(21, 22), name = "") +
  annotate(label = paste0("p=", filter(glm_res, grepl("Clust2", term)) %>% pull(p.value) %>% round(digits = 4)), 
           geom = "text", y = 5.2, x = 2.5) + theme(legend.position = "bottom")
ggsave(total_rm_plot, height = 3.8, width = 3.5, file = "total_rm_systems_clust2_boxplot_isolatesMags.pdf")    

# Split up in the other order
total_rm_plot <- ggplot(rm_counts_per_genome, 
                          aes(x = factor(Clust2), y = jitter(V1))) + 
  geom_beeswarm(method = "compactswarm", cex = 1.5,  aes(fill = factor(Clust2), 
                                                         shape = ifelse(IsMAG, "MAG", "Isolate")), alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, fill = NA) + ylab("Number of RM systems") +
  scale_x_discrete(name = "Genome group", labels = c("A", "B")) + 
  facet_grid(~ifelse(IsMAG, "MAGs", "Isolates")) + 
  scale_fill_brewer(palette = "Paired", guide = "none") + scale_shape_manual(values = c(21, 22), name = "") +
  annotate(label = paste0("Group B p=", filter(glm_res, grepl("Clust2", term)) %>% pull(p.value) %>% round(digits = 4)), 
           geom = "text", y = 5.2, x = 2) + theme(legend.position = "bottom") +
  theme(strip.background = element_rect(color = "white"))
ggsave(total_rm_plot, height = 3.8, width = 3.5, file = "total_rm_systems_clust2_boxplot_isolatesMags.pdf")    

## Test isolates alone
rm_counts_per_genome[IsMAG==F, wilcox.test(V1[Clust2==1], V1[Clust2==2])]
# p = 0.00063
rm_counts_per_genome[Clust2 ==2 & IsMAG==F, table(V1 > 1)]

rm_counts_sub <- df_data[type == "RM", length(sys_id), by = list(genome, Clust2, IsMAG, subtype)]
rm_counts_sub <- merge(rm_counts_sub, good_genomes[,list(genome, Clust2, IsMAG, Genome_Size, n_genes_called, Completeness)], all.y = T)
rm_counts_sub <- rm_counts_sub %>% dcast(genome+Clust2+IsMAG~subtype, value.var = "V1", fill = 0) %>% 
  melt(id.var = c("genome", "Clust2", "IsMAG"), variable.name = "subtype", variable.factor = F)

rm_counts_sub[is.na(value), value:=0]
rm_counts_sub <- rm_counts_sub[subtype != "NA"]

stats_rm_subtypes <- rm_counts_sub %>% 
  merge(good_genomes[,list(genome, Genome_Size)]) %>% 
  group_by(subtype) %>% 
  reframe(broom::tidy(glm(value~factor(IsMAG)+factor(Clust2)+offset(log(Genome_Size)), family = "poisson"))) %>% 
  group_by(term) %>% 
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  ungroup()
stats_rm_subtypes %>% filter(p.adj < 0.1)

rm_dist <- ggplot(rm_counts_sub[IsMAG==F], aes(x = factor(Clust2), y = factor(value))) + 
  geom_bin_2d(color = "black", breaks=seq(-0.5,5.5,1)) + facet_grid(.~subtype) +
  scale_fill_distiller(palette = "Purples", name = "Number of\nisolate genomes", direction = 1) + 
  xlab("Genome group") + 
  ylab("Copy number of RM systems") + scale_y_discrete(limits = factor(seq(0, 5))) + 
  scale_x_discrete(limits = factor(c(1,2)), labels = c("A", "B"))+
  geom_text(data = filter(stats_rm_subtypes, grepl("Clust2", term) & p.adj < 0.1), x = 1.5, y = 6, 
            aes(label = paste0("p=", format(p.adj, digits = 3))))
ggsave(rm_dist, file = "~/Desktop/globus_transfers_elGenoems/rm_subtype_dist.pdf", width = 7, height = 3.4)



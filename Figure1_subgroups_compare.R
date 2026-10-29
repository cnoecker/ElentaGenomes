library(tidyverse)
library(vegan)
library(ape)
library(ggtree)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
library(data.table)
theme_set(theme_classic2())
## Figure 1, Figure S2, Figure S3 (pangenome and phylogeography)

# Read in metadata
metadata_qc <- readRDS("final_qc_metadata.rds")
el_metadata <- metadata_qc[grepl("lenta", species)]

# Clustering and Mash distance PCoA
# Figure 1 B - PCoA Mash distances
mash_dists <- fread("all_mash_dists_good_genomes.txt")
mash_dists <- mash_dists[V1 %in% metadata_qc[,genome] & V2 %in% metadata_qc[,genome]]

el_mash_dist <- mash_dists[V1 %in% metadata_qc[grepl("lenta", species), genome] & V2 %in% metadata_qc[grepl("lenta", species), genome]]
el_dist_mat <- dcast(el_mash_dist, V1~V2, value.var = "V3") %>% column_to_rownames("V1") %>% as.matrix() %>% as.dist()
pcoa3 <- pcoa(el_dist_mat)
pcoa3_res <- pcoa3$vectors %>% data.frame() %>% rownames_to_column("genome") %>% data.table() %>% 
  merge(metadata_qc[,list(genome, Continent, Country,  IsMAG, GenomeSource.x, Clust2)], by = "genome")

continent_test <- adonis2(el_dist_mat~Continent, data = metadata_qc[grepl("lenta", species)])
pcoa_plot_mash <- ggplot(pcoa3_res, aes(x = Axis.1, y = Axis.2, fill = Continent)) + 
  geom_point(aes(shape = ifelse(IsMAG == TRUE, "MAG", "Isolate")), size = 1.8) + 
  scale_shape_manual(values = c(21,22), name = "") + 
  scale_fill_manual(values = continent_colors, guide = guide_legend(ncol = 2, override.aes = list(shape = 21))) + 
  theme_classic2() + 
  geom_text(x = 0.005, y = -0.008, label = paste0("R2 = ", round(continent_test[3][1,], 2), "\n(p=", 
                                                  continent_test[5][1,], ")"), size = 3) + ylim(-0.009, 0.013)

pcoa_plot_mash_clust <- ggplot(pcoa3_res, aes(x = Axis.1, y = Axis.2, fill = factor(Clust2))) + 
  geom_point(aes(shape = ifelse(IsMAG == TRUE, "MAG", "Isolate")), size = 1.8) + 
  scale_shape_manual(values = c(21,22), name = "") + 
  scale_fill_manual(values = brewer.pal(3, "Paired")[1:2], name = "Genome Group", guide = guide_legend(ncol = 2, override.aes = list(shape = 21))) + 
  theme_classic2()
ggsave(pcoa_plot_mash_clust + ggtitle("Mash distances"), file = "pcoa_mash_clust.pdf", width = 4.8, height = 3.2)

## Compare clustering options
clust1 <- hclust(el_dist_mat, method = "ward.D2")
clust_sub <- cutree(clust1, h = 0.06)
clust_sub2 <- cutree(clust1, h = 0.05)
clust_sub3 <- cutree(clust1, h = 0.065)

clust1 %>% plot(labels = clust_sub)
clust2 <- cutree(clust1, h = 0.08)
# Put onto PCoa
# moving forward with Clust 2 
pcoa3_res <- merge(pcoa3_res, data.table(genome = names(clust_sub), Clust = clust_sub, Clust2 = clust2, Clust3 = clust_sub3, Clust4 = clust_sub2), by = "genome")
clust_mash_pcoa <- ggplot(pcoa3_res, aes(x = Axis.1, y = Axis.2, fill = factor(Clust2))) + geom_point(aes(shape = ifelse(IsMAG == TRUE, "MAG", "Isolate")), size = 1.7) + scale_shape_manual(values = c(21,22), name = "") + 
  scale_fill_brewer(palette = "Paired", guide = guide_legend(ncol = 2, override.aes = list(shape = 21)), name = "Genome Cluster")

# Figure S3A-B
library(factoextra)
clust_hc <- hcut(el_dist_mat, k = 10)
clust_hc2 <- hcut(el_dist_mat, k = 3)
fviz_silhouette(clust_hc)
fviz_dend(clust_hc2)
fviz_dend(clust_hc)
clustplot_wss <- fviz_nbclust(as.matrix(el_dist_mat), FUN = hcut, method = "wss")
clustplot_sil <- fviz_nbclust(as.matrix(el_dist_mat), FUN = hcut, method = "silhouette")
ggsave(clustplot_wss+clustplot_sil, file = "cluster_plots_eval.pdf", width = 7, height = 3)


## Figure 1a
### Phylogenetic tree plot
# Read in tree
elenta_tree <- read.tree("RAxML_bestTree.mgs_fG2akJovke_refined.tre")

# Fix tree labels
elenta_tree$tip.label <- gsub("'", "", elenta_tree$tip.label)

# Root tree with E. timonensis as the outgroup
elenta_tree <- root(elenta_tree, "GCF_900184265.1_PRJEB20940_genomic") #timonensis
# Remove E. timonesis so that we just have E. lenta genomes
elenta_tree <- drop.tip(elenta_tree, "GCF_900184265.1_PRJEB20940_genomic")

# Set up a color scale for different continents and countries
country_tab <- metadata_qc %>% select(Continent, Country) %>% distinct() %>% 
  arrange(Continent)
continent_colors <- c(brewer.pal(5, "Set1"), "gray")
names(continent_colors) <- unique(country_tab$Continent)
num_colors <- country_tab %>% count(Continent) #Number of countries per continent
country_colors <- c(colorRampPalette(c(brewer.pal(9, "Reds")[6], brewer.pal(9, "Reds")[8]))(num_colors[1,]$n), 
                    colorRampPalette(c(brewer.pal(9, "Blues")[5], brewer.pal(9, "Blues")[9]))(num_colors[2,]$n),
                    colorRampPalette(c(brewer.pal(9, "Purples")[5], brewer.pal(9, "Purples")[9]))(num_colors[3,]$n),
                    colorRampPalette(c(brewer.pal(9, "Greens")[5], brewer.pal(9, "Greens")[9]))(num_colors[4,]$n),
                    colorRampPalette(c(brewer.pal(9, "Oranges")[3], brewer.pal(9, "Oranges")[7]))(num_colors[5,]$n), "gray30", "gray")
names(country_colors) <- country_tab$Country

metadata_qc <- metadata_qc %>% mutate(
  treeLab2=ifelse(grepl("0024265", genome), "DSM 2243", ""))

# Add a heatmap/rectangles around the tree to show the 2 subspecies
tree_viz_clust_continent <- ggtree(elenta_tree, layout = "fan") %<+% 
  select(metadata_qc, genome, species, Species2, Country, Continent, age, age_category, treeLab, treeLab2, NewGenome, Clust2) + #, Clust)] + 
  geom_hilight(node = 386, fill = colorspace::lighten(brewer.pal(3, "Paired")[2], amount = 0.25), alpha = 0.5) +
  geom_tree(aes(color = Continent), size = 0.5, layout = "fan") + 
  scale_color_manual(values = continent_colors, name = "Continent")+ 
  geom_tippoint(aes(size = ifelse(NewGenome, "New from this study", "")), color = "black", fill = brewer.pal(12, "Set3")[12], shape = 21) + 
  geom_treescale() + scale_size_manual(values = c(0.01, 2), name = "")+ 
  scale_fill_manual(values = c("gold", "lightblue", "lightgreen")) + 
  geom_tiplab(aes(label = treeLab2))
ggsave(tree_viz_clust_continent, file = "elenta_tree_clust_continent.pdf", width = 8, height = 7)


## Figure 1C and Fig S3 - ANI

el_ani <- fread("pyani_elenta/ANIm_percentage_identity.tab")
el_ani <- melt(el_ani, id.var = "V1", variable.name = "genome2", variable.factor = F)
el_ani <- el_ani[genome2 != "GCF_900184265.1_PRJEB20940_genomic" & V1 != "GCF_900184265.1_PRJEB20940_genomic"]
el_ani[,Dist:=1-value]
el_ani_dist <- dcast(el_ani, V1 ~ genome2, value.var = "Dist") %>% column_to_rownames("V1") %>% as.matrix() %>% as.dist()
el_ani_pcoa <- cmdscale(el_ani_dist) %>% data.frame() %>% rownames_to_column("genome") %>% data.table()
el_ani_pcoa <- merge(el_ani_pcoa, metadata_qc[,list(genome, Clust2, Continent, Country, IsMAG)], by = "genome", all.x = T)

adonis2(el_ani_dist~Country, data = metadata_qc[grepl("lenta", species)])
adonis2(el_ani_dist~Continent, data = metadata_qc[grepl("lenta", species)])

##PCOA following others
gene_pcoa_el <- pcoa(el_ani_dist)
gene_pcoa_res_el <- gene_pcoa_el$vectors %>% data.frame() %>% rownames_to_column("genome") %>% data.table() %>% 
  merge(metadata_qc[,list(genome, Continent, Country, species, Species2,Clust2, IsMAG)], by = "genome")
pcoa_plot_ani_el <- ggplot(gene_pcoa_res_el, aes(x = Axis.1, y = Axis.2, fill = factor(Clust2))) + 
  geom_point(aes(shape = ifelse(IsMAG == TRUE, "MAG", "Isolate")), size = 1.7) + 
  scale_shape_manual(values = c(21,22), name = "") + 
  scale_fill_manual(values = brewer.pal(3, "Paired")[1:2], name = "Genome Group", guide = guide_legend(ncol = 2, override.aes = list(shape = 21))) + 
  theme_classic2() #, name = "Country" 
ggsave(pcoa_plot_ani_el + ggtitle("Whole-genome ANI"), file = "pcoa_ani_clusts.pdf", width = 4.8, height = 3.2)

continent_test <- adonis2(el_ani_dist~Continent, data = metadata_qc[grepl("lenta", species)])
pcoa_plot_ani_el <- ggplot(gene_pcoa_res_el, aes(x = Axis.1, y = Axis.2, fill = Continent)) + 
  geom_point(aes(shape = ifelse(IsMAG == TRUE, "MAG", "Isolate")), size = 1.8) + 
  scale_shape_manual(values = c(21,22), name = "") + 
  scale_fill_manual(values = continent_colors, guide = guide_legend(ncol = 2, override.aes = list(shape = 21))) + 
  theme_classic2()  + geom_text(x = 0.005, y = -0.01, label = paste0("R2 = ", round(continent_test[3][1,], 2), "\n(p=", continent_test[5][1,], ")"), size = 3)#, name = "Country" 
ggsave(pcoa_plot_ani_el + ggtitle("Whole-genome ANI"), file = "pcoa_ani_continent.pdf", width = 5.5, height = 3.2)

el_ani_nmds <- metaMDS(el_ani_dist, k = 2, trymax = 30) %>% scores() %>% data.frame() %>% 
  rownames_to_column("genome") %>% data.table()
el_ani_nmds <- merge(el_ani_nmds, metadata[,list(genome, Clust2, Continent, Country, IsMAG)], by = "genome", all.x = T)

el_ani <- el_ani[genome2 != V1]
el_ani <- merge(el_ani, metadata_qc[,list(genome, Clust2, Continent, Country, IsMAG)], by.x = "V1", by.y = "genome", all.x = T)
el_ani <- merge(el_ani, metadata_qc[,list(genome, Clust2, Continent, Country, IsMAG)], by.x = "genome2", by.y = "genome", all.x = T)

el_ani[,Comparison:=case_when(
  Clust2.x == 1 & Clust2.y == 1 ~ "Within\nGroup 1", 
  Clust2.x == 2 & Clust2.y == 2 ~ "Within\nGroup 2", 
  TRUE ~ "Between\nGroups"
)]
ani_compare_plot <- ggplot(el_ani[genome2 < V1], aes(y = value, x = Comparison, fill = Comparison)) + geom_violin(draw_quantiles = c(0.5)) +
  ylab("Average nucleotide identity") + theme(axis.text = element_text(size = 11))+ 
  scale_fill_manual(values = c("gray", brewer.pal(5, "Blues")[c(2,5)]), guide = "none") +
  xlab("") + stat_compare_means(comparisons = list(c(1,2), c(1,3), c(2,3)), symnum.args = list(cutpoints = c(0, 0.001, Inf), symbols = c("***", "ns")))
#+
ggsave(ani_compare_plot, file = "ani_dist_plot.pdf", width = 3.6, height = 4.2)

el_cov <- fread("pyani_elenta/ANIm_alignment_coverage.tab")
el_cov <- melt(el_cov, id.var = "V1", variable.name = "genome2", variable.factor = F)
el_cov <- el_cov[genome2 != "GCF_900184265.1_PRJEB20940_genomic" & V1 != "GCF_900184265.1_PRJEB20940_genomic"]
setnames(el_cov, "value", "cov")
el_ani <- merge(el_ani, el_cov, by = c("V1", "genome2"))

## Figure S3D - Heatmap
library(ComplexHeatmap)
clust_annot <- metadata[,list(genome, Group = as.character(Clust2))] %>% column_to_rownames("genome")
clust_cols <- brewer.pal(3, "Paired")[1:2]
names(clust_cols) <- c("1", "2")
pdf("ani_heatmap_clust.pdf", width = 8, height = 7)
el_ani %>% dcast(V1 ~ genome2, value.var = "value", fill = 1) %>% column_to_rownames("V1") %>% as.matrix() %>% 
  Heatmap(row_labels = rep("", 284), column_labels = rep("", 284), 
          heatmap_legend_param = list(title = "ANI"), column_title = "Genomes", row_title = "Genomes", 
          col = circlize::colorRamp2(breaks = c(0.97, 0.987, 1), colors = brewer.pal(11, "RdGy")[c(7,6,1)]),
          top_annotation = HeatmapAnnotation(df = clust_annot, col = list(Group = clust_cols[1:2])))
dev.off()

## Figure S3 F, I - Gene family pres/abs PCoA
# % in accessory genome per genome
eg_el_annots <- merge(eg_el_annots, el_metadata[,list(genome, IsMAG, Clust2)], by="genome", all.x = T)
eg_el_annots[,1-sum(SoftCore)/length(SoftCore), by=genome][,summary(V1)]
eg_el_annots[IsMAG == F,1-sum(SoftCore)/length(SoftCore), by=genome][,summary(V1)]
eggnog_isolates <- eg_el_annots[IsMAG==F]
isolate_counts <- eggnog_isolates[,length(unique(genome)), by= eggNOG_OGs]
eggnog_isolates[,CoreIsolate:=ifelse(eggNOG_OGs %in% isolate_counts[V1 == 168, eggNOG_OGs], 1, 0)]
eggnog_isolates[,1-sum(CoreIsolate)/length(CoreIsolate), by=genome][,summary(V1)]

elenta_wide <- dcast(eg_el_annots, genome~ eggNOG_OGs, value.var = "Present", fun.aggregate = sum)
genome_dist_elenta <- vegdist(as.matrix(column_to_rownames(elenta_wide, "genome")), metric = "manhattan")
gene_pcoa_el <- pcoa(genome_dist_elenta)
#Continent and clust comparison
gene_pcoa_res_el <- gene_pcoa_el$vectors %>% data.frame() %>% rownames_to_column("genome") %>% data.table() %>% merge(metadata_qc[,list(genome, Continent, Country, species, Species2, Clust, Clust2, IsMAG)], by = "genome")

continent_test <- adonis2(genome_dist_elenta~Continent, data = metadata_qc[grepl("lenta", species)])
pcoa_plot_genes_el <- ggplot(gene_pcoa_res_el, aes(x = Axis.1, y = Axis.2, fill = Continent)) + 
  geom_point(aes(shape = ifelse(IsMAG == TRUE, "MAG", "Isolate")), size = 1.8) + 
  scale_shape_manual(values = c(21,22), name = "") + 
  scale_fill_manual(values = continent_colors, guide = guide_legend(ncol = 2, override.aes = list(shape = 21))) + 
  theme_classic2()  + geom_text(x = 0.06, y = -0.06, label = paste0("R2 = ", round(continent_test[3][1,], 2), "\n(p=", continent_test[5][1,], ")"), size = 3)#, name = "Country" #, name = "Country"
ggsave(pcoa_plot_genes_el + ggtitle("Gene Family Content"), file = "pcoa_genefamily_continent.pdf", width = 5.5, height = 3.2)

pcoa_plot_genes_el_clust <- ggplot(gene_pcoa_res_el, aes(x = Axis.1, y = Axis.2, fill = factor(Clust2))) + 
  geom_point(aes(shape = ifelse(IsMAG == TRUE, "MAG", "Isolate")), size = 1.7) + 
  scale_shape_manual(values = c(21,22), name = "") + 
  scale_fill_manual(values = brewer.pal(3, "Paired")[1:2], name = "Genome Cluster", guide = guide_legend(ncol = 2, override.aes = list(shape = 21))) + 
  theme_classic2() #, name = "Country" 
ggsave(pcoa_plot_genes_el_clust + ggtitle("Gene Family Content"), file = "pcoa_genefamily_clust.pdf", width = 4.8, height = 3.2)

## Figure 1D
el_metadata[,fisher.test(Clust2, Continent)]

continent_dist <- ggplot(el_metadata, aes(x = Clust2, fill = Continent)) + geom_bar(position = "fill", color = "black") + 
  scale_fill_manual(values = continent_colors) + theme_classic2() + xlab("Genome group") + 
  ylab("Fraction of genomes")
ggsave(continent_dist, file = "continent_clust_compare.pdf", width = 4.5, height = 5)

## Genome source compare - Figure S3J
el_metadata[,GenomeSource2:=factor(case_when(
  GenomeSource.x %in% c("Other", "Harris Wang lab") ~ "Other",
  GenomeSource.x == "ImmunoMicrobiome" ~ "This study",
  TRUE ~ GenomeSource.x
), levels = c("This study", "NCBI", "UHGG", "Other"))]
genome_source_dist <- ggplot(el_metadata, aes(x = Clust2, fill = GenomeSource2)) + geom_bar(position = "fill", color = "black") + 
  scale_fill_brewer(palette = "Dark2", name = "Genome source") + theme_classic2() + xlab("Genome group") + 
  ylab("Fraction of genomes")
ggsave(genome_source_dist, file = "genomesource_clust_compare.pdf", width = 4.5, height = 5)
el_metadata[,fisher.test(GenomeSource2, Clust2)]

## Figure S2A
# Make a separate summary of isolate/mag sources
genome_sources <- ggplot(metadata_qc, aes(x = ifelse(NewGenome, "From this study", "Other"), fill = Genome_type)) +
  geom_bar(color = "black") + theme_classic2() + 
  scale_fill_brewer(palette = "Set2", name = "") +
  xlab("") + ylab("Number of genomes") + theme(legend.position ="bottom")
ggsave(genome_sources, file = "elenta_genome_set_summary.pdf", width = 2.5, height = 4)

# S2B-C - switching to CheckM2 to match the supp table
mag_qc <- el_metadata %>% filter(Genome_type=="MAG") %>% select(Completeness, Contamination) %>% 
  pivot_longer(cols = 1:2) %>% 
  ggplot(aes(x = name, y = value))+ geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.4, height = 0, size = 0.5) + facet_wrap(~name, scales = "free") + ylab("Percent") + xlab("") +
  theme_classic2() + theme(axis.text.x = element_blank())
ggsave(mag_qc, file = "mag_checkm2_completeness_contam.pdf", width = 3.5, height = 4)

# S2D-E Pangenome sampling
eg_el_annots <- readRDS("elenta_eggnog.rds")

### Pan-genome sampling
genomes <- unique(eg_el_annots[,genome])
#eg_counts_el <- eg_el_annots[,length(`#query`), by = list(eggNOG_OGs, seed_ortholog, genome)]
eg_counts_el <- eg_el_annots[,length(`#query`), by = list(eggNOG_OGs,genome)]
eg_counts_el[,Present:=1]
OGs <- eg_counts_el 

Pangenome<-data.table(Ngenomes=numeric(), Iteration=numeric(), New=numeric(), Pangenome=numeric(), Core=numeric(),Accessory=numeric())

# Ngis the number of genomes to plot, and i is the sampling
Pangenome <- rbindlist(lapply(2:length(genomes), function(Ng){
  message("Ngenome=", Ng)
  tab_sub <- data.table(Ngenome = rep(Ng), Iteration = 1:length(genomes),
                        New = rep(NA), Pangenome = rep(NA), Core = rep(NA), SoftCore = rep(NA), Accessory = rep(NA))
  tab_sub <- rbindlist(lapply(1:length(genomes), function(i){
    sgenomes <- sample(genomes, Ng, replace=F)
    
    new_genes <- OGs[genome %in% sgenomes]
    new_genes[,New:=ifelse(genome == sgenomes[Ng], TRUE, FALSE)] # If it is present in the last sampled genome (vs. other genomes)
    new_obs <- new_genes[,sum(Present), by = list(New, eggNOG_OGs)] %>% dcast(eggNOG_OGs~New, value.var = "V1", fill = 0)
    new_obs[,NewObs:=ifelse(`TRUE`==1 & `FALSE`==0, T, F)]
    new_obs[,TotGenomes:=`FALSE`+`TRUE`]
    num_new <- new_obs[,sum(NewObs)]
    foo <- data.table(Ngenome = Ng, Iteration = i, New = num_new,
                      Pangenome = nrow(new_obs),
                      Core = nrow(new_obs[TotGenomes==Ng]),
                      SoftCore = nrow(new_obs[TotGenomes >= 0.95*Ng]),
                      Accessory = nrow(new_obs[TotGenomes < Ng]), 
                      SoftAccessory = nrow(new_obs[TotGenomes < 0.95*Ng]))
    return(foo)
  }))
  return(tab_sub)
}), fill = T)

saveRDS(Pangenome, file = "pangenome_sampling_results_softcore_eggNOGOGs.RDS")

# Plot curves
pangenome_curve_plot <- Pangenome %>%
  melt(id.vars = c("Ngenome", "Iteration"), value.name = "Count", variable.name = "Metric") %>% 
  filter(Metric %in% c("SoftCore", "Accessory", "New", "Pangenome")) %>% 
  #gather(-Ngenomes, -Iteration, key=Metric, value=Count) %>%
  ggplot(aes(x=Ngenome, y=Count, fill=factor(Metric, levels = c("Pangenome", "Accessory", "SoftCore", "New")))) +
  stat_summary(geom="ribbon", fun.data=mean_sd, alpha=0.5) +
  theme_classic2() +
  xlab("Number of genomes analyzed") +
  ylab("Number of gene families") +
  xlim(0, 284) + scale_fill_brewer(palette = "Dark2", name = "", labels = c("Pangenome", "Accessory", "Core (>95%)", "New"))

ggsave(pangenome_curve_plot, file = "pangenome_curve_Elenta_all_eggNOGOGs.pdf", width = 4, height = 2.5)

Pangenomestats<-
  Pangenome %>%
  select(Ngenome, New, Pangenome, Core, Accessory) %>%
  melt(id.var = "Ngenome", variable.name = "Metric", value.name = "Count") %>%
  group_by(Ngenome, Metric) %>%
  summarize_at("Count", list(mean = mean,median = median,min = min,max = max)) %>%
  arrange(desc(Ngenome))

### Isolates only
genomes <- unique(eg_el_annots[genome %in% metadata_qc[IsMAG == F, genome],genome])
eg_counts_el <- eg_el_annots[,length(`#query`), by = list(eggNOG_OGs, genome)]
eg_counts_el <- eg_counts_el[genome %in% genomes]
eg_counts_el[,Present:=1]
OGs <- eg_counts_el 

Pangenome<-data.table(Ngenomes=numeric(), Iteration=numeric(), New=numeric(), Pangenome=numeric(), Core=numeric(),Accessory=numeric())

# Ngis the number of genomes to plot, and i is the sampling
Pangenome <- rbindlist(lapply(2:length(genomes), function(Ng){
  message("Ngenome=", Ng)
  tab_sub <- data.table(Ngenome = rep(Ng), Iteration = 1:length(genomes),
                        New = rep(NA), Pangenome = rep(NA), Core = rep(NA), SoftCore = rep(NA), Accessory = rep(NA))
  tab_sub <- rbindlist(lapply(1:length(genomes), function(i){
    sgenomes <- sample(genomes, Ng, replace=F)
    new_genes <- OGs[genome %in% sgenomes]
    new_genes[,New:=ifelse(genome == sgenomes[Ng], TRUE, FALSE)] # If it is present in the last sampled genome (vs. other genomes)
    new_obs <- new_genes[,sum(Present), by = list(New, eggNOG_OGs)] %>% dcast(eggNOG_OGs~New, value.var = "V1", fill = 0)
    new_obs[,NewObs:=ifelse(`TRUE`==1 & `FALSE`==0, T, F)]
    new_obs[,TotGenomes:=`FALSE`+`TRUE`]
    num_new <- new_obs[,sum(NewObs)]
    foo <- data.table(Ngenome = Ng, Iteration = i, New = num_new,
                      Pangenome = nrow(new_obs),
                      Core = nrow(new_obs[TotGenomes==Ng]),
                      Accessory = nrow(new_obs[TotGenomes < Ng]))
    return(foo)
  }))
  return(tab_sub)
}), fill = T)

saveRDS(Pangenome, file = "pangenome_sampling_results_softcore_isolatesOnly_eggNOGOGs.RDS")

## Isolates only version
Pangenome <- readRDS("pangenome_sampling_results_softcore_isolatesOnly_eggNOGOGs.RDS")

pangenome_curve_plot <- Pangenome %>%
  melt(id.vars = c("Ngenome", "Iteration"), value.name = "Count", variable.name = "Metric") %>% 
  filter(Metric %in% c("Core", "Accessory", "New", "Pangenome")) %>% 
  #gather(-Ngenomes, -Iteration, key=Metric, value=Count) %>%
  ggplot(aes(x=Ngenome, y=Count, fill=factor(Metric, levels = c("Pangenome", "Accessory", "Core", "New")))) +
  stat_summary(geom="ribbon", fun.data=mean_sd, alpha=0.5) +
  theme_classic2() +
  xlab("Number of isolate genomes analyzed") +
  ylab("Number of gene families") +
  xlim(0, 168) + scale_fill_brewer(palette = "Dark2", name = "", labels = c("Pangenome", "Accessory", "Core", "New"))
ggsave(pangenome_curve_plot, file = "pangenome_curve_Elenta_isolates_eggNOGOGs.pdf", width = 4, height = 2.5)

# Heaps law
library(micropan)
eg_el_annots[,EGClean:=gsub("[@,\\|]", "_", eggNOG_OGs)]
panmat <- dcast(eg_el_annots, genome~EGClean, value.var = "Present", fun.aggregate = sum, fill = 0) %>% 
  column_to_rownames("genome") %>% as.matrix()
heaps_mod1 <- heaps(pan.matrix = panmat, n.perm = 500)
## Isolates only
panmat_isolates <- dcast(eg_el_annots[genome %in% metadata_qc[IsMAG == F, genome]], genome~EGClean, value.var = "Present", fun.aggregate = sum, fill = 0) %>% 
  column_to_rownames("genome") %>% as.matrix()
heaps_mod2 <- heaps(pan.matrix = panmat_isolates, n.perm = 500)

## Repeat with Panaroo
panaroo_pres_abs <- fread("gene_presence_absence.csv")
panaroo_pres_abs <- melt(panaroo_pres_abs, id.var = names(panaroo_pres_abs)[1:3], variable.factor = F)
panaroo_pres_abs[,Present:=ifelse(value == "", 0, 1)]
setnames(panaroo_pres_abs, "variable", "genome")
panaroo_prev <- panaroo_pres_abs[,sum(Present), by = Gene]
panaroo_prev[,table(V1/284 > 0.95)]
core_panaroo <- panaroo_prev[V1/284 > 0.95, Gene]
panaroo_pres_abs[,SoftCore:=ifelse(Gene %in% core_panaroo, 1, 0)]
panaroo_pres_abs[Present==1,1-sum(SoftCore)/length(SoftCore), by=genome][,summary(V1)]

# ## Panaroo sampling
genomes <- unique(el_metadata[,genome])
OGs <- panaroo_pres_abs #%>% dcast(seed_ortholog~Genome, value.var = "Present", fill = 0) #%>% as.data.frame() %>% column_to_rownames("seed_ortholog")
#OGs<-OGs[rowSums(OGs)>0,]
OGs <- OGs[Present==1]
Pangenome<-data.table(Ngenomes=numeric(), Iteration=numeric(), New=numeric(), Pangenome=numeric(), Core=numeric(),Accessory=numeric())

# Ngis the number of genomes to plot, and i is the sampling
Pangenome <- rbindlist(lapply(2:length(genomes), function(Ng){
  message("Ngenome=", Ng)
  tab_sub <- data.table(Ngenome = rep(Ng), Iteration = 1:length(genomes),
                        New = rep(NA), Pangenome = rep(NA), Core = rep(NA), SoftCore = rep(NA), Accessory = rep(NA))
  tab_sub <- rbindlist(lapply(1:length(genomes), function(i){
    sgenomes <- sample(genomes, Ng, replace=F)

    new_genes <- OGs[genome %in% sgenomes]
    new_genes[,New:=ifelse(genome == sgenomes[Ng], TRUE, FALSE)] # If it is present in the last sampled genome (vs. other genomes)
    new_obs <- new_genes[,sum(Present), by = list(New, Gene)] %>% dcast(Gene~New, value.var = "V1", fill = 0)
    new_obs[,NewObs:=ifelse(`TRUE`==1 & `FALSE`==0, T, F)]
    new_obs[,TotGenomes:=`FALSE`+`TRUE`]
    num_new <- new_obs[,sum(NewObs)]
    foo <- data.table(Ngenome = Ng, Iteration = i, New = num_new,
                      Pangenome = nrow(new_obs),
                      Core = nrow(new_obs[TotGenomes==Ng]),
                      SoftCore = nrow(new_obs[TotGenomes >= 0.95*Ng]),
                      Accessory = nrow(new_obs[TotGenomes < Ng]),
                      SoftAccessory = nrow(new_obs[TotGenomes < 0.95*Ng]))
    return(foo)
  }))
  return(tab_sub)
}), fill = T)

saveRDS(Pangenome, file = "pangenome_sampling_results_softcore_panaroo.RDS")

## Panaroo sampling
genomes <- unique(el_metadata[IsMAG==F,genome])
OGs <- panaroo_pres_abs 
OGs <- OGs[Present==1 & genome %in% genomes]
Pangenome<-data.table(Ngenomes=numeric(), Iteration=numeric(), New=numeric(), Pangenome=numeric(), Core=numeric(),Accessory=numeric())

# Ngis the number of genomes to plot, and i is the sampling
Pangenome <- rbindlist(lapply(2:length(genomes), function(Ng){
  message("Ngenome=", Ng)
  tab_sub <- data.table(Ngenome = rep(Ng), Iteration = 1:length(genomes),
                        New = rep(NA), Pangenome = rep(NA), Core = rep(NA), SoftCore = rep(NA), Accessory = rep(NA))
  tab_sub <- rbindlist(lapply(1:length(genomes), function(i){
    sgenomes <- sample(genomes, Ng, replace=F)
    
    new_genes <- OGs[genome %in% sgenomes]
    new_genes[,New:=ifelse(genome == sgenomes[Ng], TRUE, FALSE)] # If it is present in the last sampled genome (vs. other genomes)
    new_obs <- new_genes[,sum(Present), by = list(New, Gene)] %>% dcast(Gene~New, value.var = "V1", fill = 0)
    new_obs[,NewObs:=ifelse(`TRUE`==1 & `FALSE`==0, T, F)]
    new_obs[,TotGenomes:=`FALSE`+`TRUE`]
    num_new <- new_obs[,sum(NewObs)]
    foo <- data.table(Ngenome = Ng, Iteration = i, New = num_new,
                      Pangenome = nrow(new_obs),
                      Core = nrow(new_obs[TotGenomes==Ng]),
                      SoftCore = nrow(new_obs[TotGenomes >= 0.95*Ng]),
                      Accessory = nrow(new_obs[TotGenomes < Ng]), 
                      SoftAccessory = nrow(new_obs[TotGenomes < 0.95*Ng]))
    return(foo)
  }))
  return(tab_sub)
}), fill = T)

saveRDS(Pangenome, file = "pangenome_sampling_results_softcore_panaroo_isolatesOnly.RDS")

panaroo_samp <- readRDS("pangenome_sampling_results_softcore_panaroo.RDS")
## All resampling curves
pangenome_curve_plot <- panaroo_samp %>%
  melt(id.vars = c("Ngenome", "Iteration"), value.name = "Count", variable.name = "Metric") %>% 
  filter(Metric %in% c("SoftCore", "SoftAccessory", "New", "Pangenome")) %>% 
  #gather(-Ngenomes, -Iteration, key=Metric, value=Count) %>%
  ggplot(aes(x=Ngenome, y=Count, fill=factor(Metric, levels = c("Pangenome", "SoftAccessory", "SoftCore", "New")))) +
  stat_summary(geom="ribbon", fun.data=mean_sd, alpha=0.5) +
  theme_classic2() +
  xlab("Number of isolate genomes analyzed") +
  ylab("Number of gene families") +
  scale_fill_brewer(palette = "Dark2", name = "", labels = c("Pangenome", "Accessory", "Core", "New"))
ggsave(pangenome_curve_plot, file = "pangenome_curve_Elenta_panaroo.pdf", width = 4, height = 2.5)

panaroo_samp_isolates <- readRDS("pangenome_sampling_results_softcore_panaroo_isolatesOnly.RDS")

pangenome_curve_plot <- panaroo_samp_isolates %>%
  melt(id.vars = c("Ngenome", "Iteration"), value.name = "Count", variable.name = "Metric") %>% 
  filter(Metric %in% c("Core", "Accessory", "New", "Pangenome")) %>% 
  #gather(-Ngenomes, -Iteration, key=Metric, value=Count) %>%
  ggplot(aes(x=Ngenome, y=Count, fill=factor(Metric, levels = c("Pangenome", "Accessory", "Core", "New")))) +
  stat_summary(geom="ribbon", fun.data=mean_sd, alpha=0.5) +
  theme_classic2() +
  xlab("Number of isolate genomes analyzed") +
  ylab("Number of gene families") +
  scale_fill_brewer(palette = "Dark2", name = "", labels = c("Pangenome", "Accessory", "Core", "New"))
ggsave(pangenome_curve_plot, file = "pangenome_curve_Elenta_panaroo_isolatesOnly.pdf", width = 4, height = 2.5)

# Heaps law
panmat <- dcast(panaroo_pres_abs, genome~Gene, value.var = "Present", fun.aggregate = sum, fill = 0) %>% 
  column_to_rownames("genome") %>% as.matrix()
heaps_mod1 <- heaps(pan.matrix = panmat, n.perm = 500) #alpha = 0.928
## Isolates only
panmat_isolates <- dcast(panaroo_pres_abs[genome %in% el_metadata[IsMAG==F, genome] & Present==1], genome~Gene, value.var = "Present", fun.aggregate = sum, fill = 0) %>% 
  column_to_rownames("genome") %>% as.matrix()
heaps_mod2 <- heaps(pan.matrix = panmat_isolates, n.perm = 500) #alpha = 0.91


# Figure S2F - pangenome categories
top_cat <- eg_el_annots[,length(unique(eggNOG_OGs)), by=COG_category][order(V1, decreasing = T)][1:10, COG_category]
eg_el_annots[,COG_category2:=ifelse(COG_category %in% top_cat, COG_category, "Multiple/Other")]

cat_names <- data.table(COG_category = c(eg_el_annots[nchar(COG_category)==1, sort(unique(COG_category))], "W", "Z"), 
                        CategoryName =c(
                          "Not assigned",
                          #"Chromatin structure and dynamics",
                          "Energy production and conversion",
                          "Cell cycle control, cell division, chromosome partitioning",
                          "Amino acid transport and metabolism",
                          "Nucleotide transport and metabolism",
                          "Carbohydrate transport and metabolism",
                          "Coenzyme transport and metabolism",
                          "Lipid transport and metabolism",
                          "Translation, ribosomal structure and biogenesis",
                          "Transcription",
                          "Replication, recombination and repair",
                          "Cell wall/membrane/envelope biogenesis",
                          "Cell motility",
                          "Posttranslational modification, protein turnover, chaperones",
                          "Inorganic ion transport and metabolism",
                          "Secondary metabolites biosynthesis, transport and catabolism",
                          "Function unknown",
                          "Signal transduction mechanisms",
                          "Intracellular trafficking, secretion, and vesicular transport",
                          "Defense mechanisms",
                          "Extracellular structures",
                          "Cytoskeleton"))

# Fresh data setup
eg_el_annots <- readRDS("elenta_eggnog.rds")
eg_el_annots <- merge(eg_el_annots, metadata_qc[,list(genome, IsMAG, species, Species2, Country, Continent, birthCountryMerged, birthContinentMerged, 
                                                        age, gender, AnimalHost, BodySite, AnimalHost2, Clust, Clust2)])
clust_sets <- unique(eg_el_annots[,list(eggNOG_OGs, COG_category, Description)])
clust_sets[,Clust:=sapply(eggNOG_OGs, function(x){
  return(sort(unique(eg_el_annots[eggNOG_OGs==x, Clust2])))
})]
clust_sets[,Clust1:=sapply(Clust, function(x){ if (1 %in% x) return (TRUE) else return(FALSE)})]
clust_sets[,Clust2:=sapply(Clust, function(x){ if (2 %in% x) return (TRUE) else return(FALSE)})]
gene_prev <- eg_el_annots[,list(TotNum = length(unique(genome)), 
                                 TotClust1 = length(unique(genome[Clust2 == 1])), 
                                 TotClust2 = length(unique(genome[Clust2 == 2]))), by = eggNOG_OGs]
clust_sets <- merge(clust_sets, gene_prev, by = "eggNOG_OGs")
clust_sets[,Clust1_2_ratio:=TotClust2/TotClust1]
clust_sets[,LogRatio:=log2(Clust1_2_ratio) ]
clust_sets[,Taxon:=gsub(".*\\|", "", eggNOG_OGs)]
clust_sets[,SoftCore:=ifelse(TotNum > eg_el_annots[,length(unique(genome))]*0.95, 1, 0)]

clust_sets[,COG_category2:=ifelse(COG_category %in% top_cat, COG_category, "Multiple/Other")]
clust_sets <- merge(clust_sets, cat_names, by.x = "COG_category2", by.y = "COG_category", all.x = T)

clust_sets[COG_category2 %in% c("S", "-"), COG_category2:="S-"]
clust_sets[COG_category2 =="S-", CategoryName:="Unknown function/Not assigned"]
clust_sets[COG_category2=="Multiple/Other", CategoryName:=""]
cat_order <- clust_sets[,length(eggNOG_OGs), by=paste0(COG_category2, " ", CategoryName)][order(V1, decreasing = T), paste0]
category_counts <- ggplot(clust_sets, aes(x = paste0(COG_category2, " ", CategoryName), fill = ifelse(SoftCore == 1, "Core", "Accessory"))) + 
  geom_bar(color = "black")+ coord_flip() + theme_classic2() + scale_x_discrete(limits = rev(cat_order), name = "") +
  ylab("Number of gene families") + scale_fill_manual(values = brewer.pal(4, "Purples")[c(2,4)], name = "") + 
  theme(axis.text = element_text(size = 12), legend.text = element_text(size = 12))
ggsave(category_counts, file = "category_counts_core_accessory.pdf", width = 7.5, height = 4.5)

# eggnog_isolates <- eg_el_annots[IsMAG==F]
# isolate_counts <- eggnog_isolates[,length(unique(genome)), by= eggNOG_OGs]

clust_sets_isolates <- clust_sets[eggNOG_OGs %in% eggnog_isolates[,eggNOG_OGs]]
clust_sets_isolates[,CoreIsolate:=ifelse(eggNOG_OGs %in% isolate_counts[V1 == 168, eggNOG_OGs], 1, 0)]
category_counts_isolates <- ggplot(clust_sets_isolates, aes(x = paste0(COG_category2, " ", CategoryName), fill = ifelse(CoreIsolate == 1, "Core", "Accessory"))) + 
  geom_bar(color = "black")+ coord_flip() + theme_classic2() + scale_x_discrete(limits = rev(cat_order), name = "") +
  ylab("Number of gene families") + scale_fill_manual(values = brewer.pal(4, "Purples")[c(2,4)], name = "") + 
  theme(axis.text = element_text(size = 12), legend.text = element_text(size = 12))
ggsave(category_counts_isolates, file = "category_counts_isolates.pdf", width = 7.5, height = 4.5)

combined_plot <- category_counts + ggtitle("All genomes") + 
  category_counts_isolates + theme(axis.text.y = element_blank()) + ggtitle("Isolate genomes") +
  plot_layout(guides = "collect") 
ggsave(combined_plot, file = "category_counts_pangenome_combined.pdf", width = 9.5, height = 4.8)

# Figure S2G
num_enzymes <- eg_el_annots[grepl("COG0243", eggNOG_OGs)|grepl("COG1053", eggNOG_OGs), length(unique(`#query`)), 
                             by = list(genome, IsMAG, TopLevel_eggNOG)]
num_enzyme_plot <- ggplot(num_enzymes, aes(y =V1, x = ifelse(IsMAG, "MAGs", "Isolates"))) + facet_wrap(~TopLevel_eggNOG) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.25, alpha = 0.7, shape = 21, aes(fill = ifelse(IsMAG, "MAG", "Isolate"))) + ylim(0, 75) + theme_classic() +
  ylab("Number of gene copies") + xlab("") + scale_fill_brewer(palette = "Set2", guide = "none")
ggsave(num_enzyme_plot, file = "respenz_copyNum_dist.pdf", width = 5, height = 3.8)

## Figure S2H
library(panstripe)
library(ape)
# Load files
pa <- read_rtab("gene_presence_absence.Rtab")
#Filter to only isolates
pa <- pa[row.names(pa) %in% el_metadata[IsMAG==F, genome],]
tree1 <- read.tree("RAxML_bestTree.mgs_fG2akJovke_refined.tre")
tree1$tip.label <- gsub("'", "", tree1$tip.label)
# Root tree with E. timonensis as the outgroup
tree1 <- root(tree1, "GCF_900184265.1_PRJEB20940_genomic") #timonensis
tree1 <- drop.tip(tree1, "GCF_900184265.1_PRJEB20940_genomic")
tree1 <- drop.tip(tree1, el_metadata[IsMAG==T, genome])

## foo <- Biostrings::readDNAStringSet("mgs_fG2akJovke_concatenated.aln")
## width(foo) # 146058
tree1$edge.length <- tree1$edge.length*146058
# Run panstripe
fit <- panstripe(pa, tree1)
fit$summary

plot_pangenome_params(fit)
plot_pangenome_cumulative(fit)
plot_residuals(fit)
panstripe_result_tree <- plot_gain_loss(fit, tip_label = F) + 
  scale_color_gradient2(low = brewer.pal(9, "RdGy")[9], mid = brewer.pal(9, "RdGy")[4], high = brewer.pal(9, "RdGy")[1], midpoint = 300) + 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
ggsave(panstripe_result_tree, file = "panstripe_gainloss_tree.pdf", width = 9.5, height = 10.5)

panstripe_result_tree <- plot_gain_loss(fit, tip_label = T) + 
  scale_color_gradient2(low = brewer.pal(9, "RdGy")[9], mid = brewer.pal(9, "RdGy")[4], high = brewer.pal(9, "RdGy")[1], midpoint = 300) + 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

clust_assign <- el_metadata[genome %in% rownames(pa),Clust2]
names(clust_assign) <- el_metadata[genome %in% rownames(pa),genome]

panaroo_tsne <- plot_tsne(pa, category = clust_assign) + scale_color_manual(values = brewer.pal(5, "Blues")[c(2,5)], name = "Group")
ggsave(panaroo_tsne, file = "panaroo_presAbs_tSne_clust2.pdf", width = 4, height = 3.7)

## Panaroo tree - Figure S3C
plot_elenta_tree <- function(tree1){
  # Fix tree labels
  tree1$tip.label <- gsub("'", "", tree1$tip.label)
  
  if("GCF_900184265.1_PRJEB20940_genomic" %in% tree1$tip.label){
    # Root tree with E. timonensis as the outgroup
    tree1 <- root(tree1, "GCF_900184265.1_PRJEB20940_genomic") #timonensis
    # Remove E. timonesis so that we just have E. lenta genomes
    tree1 <- drop.tip(tree1, "GCF_900184265.1_PRJEB20940_genomic")
  }
  
  # Draw the tree
  elenta_tree_viz <- ggtree(tree1, layout = "fan") %<+% 
    select(el_metadata, genome, species, Species2, Country, Continent, age, age_category, treeLab, treeLab2, NewGenome, Clust2) + #, Clust)] + 
    geom_tree(aes(color = Continent), size = 0.5, layout = "fan") + 
    scale_color_manual(values = continent_colors, name = "Continent")+ 
    geom_tippoint(aes(size = ifelse(NewGenome, "New from this study", ""), shape = factor(Clust2), fill = factor(Clust2)), color = "black") + 
    geom_treescale() + scale_size_manual(values = c(0.7, 2), name = "")+ 
    scale_shape_manual(values = c(21, 22), name = "Genome group") +
    scale_fill_manual(values = c(brewer.pal(3, "Blues")[c(1,3)]), name = "Genome group") +
    geom_tiplab(aes(label = treeLab2))
  
  return(elenta_tree_viz)
}
## Panaroo tree with boostrap values
raxml_boot_tree <- read.tree("RAxML_bipartitions.raxml_boot.tre")
phylop_root <- "GCA_958428755.1_ERR9578289_bin.8_MetaWRAP_v1.3_MAG_genomic"
raxml_boot_tree <- root(raxml_boot_tree, phylop_root)
panaroo_tree_viz <- plot_elenta_tree(raxml_boot_tree) + geom_nodelab(size = 2.5)
ggsave(panaroo_tree_viz, file = "panaroo_raxml_tree_boot.pdf", width = 10, height = 8)

## Figure S3L -Birth continent PCoA
mash_dists_sub <- mash_dists[V1 %in% metadata_qc[grepl("lenta", species) & GenomeSource.x=="ImmunoMicrobiome" & !is.na(birthCountry), genome] & 
                               V2 %in% metadata_qc[grepl("lenta", species) & GenomeSource.x=="ImmunoMicrobiome" & !is.na(birthCountry), genome]]
el_dist_mat_sub <- dcast(mash_dists_sub, V1~V2, value.var = "V3") %>% column_to_rownames("V1") %>% as.matrix() %>% as.dist()
pcoa4 <- pcoa(el_dist_mat_sub)
pcoa4_res <- pcoa4$vectors %>% data.frame() %>% rownames_to_column("genome") %>% data.table() %>% 
  merge(metadata_qc[,list(genome, Continent, Country,  IsMAG, GenomeSource.x, Clust2, birthCountry, birthContinent)], by = "genome")
pcoa_plot_mash_sub <- ggplot(pcoa4_res, aes(x = Axis.1, y = Axis.2, fill = birthContinent)) + 
  geom_point(aes(shape = ifelse(IsMAG == TRUE, "MAG", "Isolate")), size = 2.1) + 
  scale_shape_manual(values = c(21,22), name = "") + 
  scale_fill_manual(values = continent_colors, guide = guide_legend(ncol = 2, override.aes = list(shape = 21)), name = "Continent of birth") + 
  theme_classic2()
ggsave(pcoa_plot_mash_sub, file = "pcoa_birthContinent_immunomicrobiomeOnly.pdf", width = 5.5, height = 3.2)
adonis2(el_dist_mat_sub~birthContinent, data = metadata_qc[grepl("lenta", species) & GenomeSource.x=="ImmunoMicrobiome" & !is.na(birthCountry)])
adonis2(el_dist_mat_sub~birthCountry, data = metadata_qc[grepl("lenta", species) & GenomeSource.x=="ImmunoMicrobiome" & !is.na(birthCountry)])

# Figure S3K - Slackia geography
mash_dists <- fread("all_mash_dists_good_genomes.txt")
mash_dists <- mash_dists[V1 %in% metadata_qc[,genome] & V2 %in% metadata_qc[,genome]]
si_mash_dist <- mash_dists[V1 %in% metadata_qc[grepl("isoflavoni", species), genome] & V2 %in% metadata_qc[grepl("isoflavoni", species), genome]]
si_dist_mat <- dcast(si_mash_dist, V1~V2, value.var = "V3") %>% column_to_rownames("V1") %>% as.matrix() %>% as.dist()
adonis2(si_dist_mat~Country, data = metadata_qc[grepl("isoflavoni", species)])
adonis2(si_dist_mat~Continent, data = metadata_qc[grepl("isoflavoni", species)])
pcoa4 <- pcoa(si_dist_mat)
pcoa4_res <- pcoa4$vectors %>% data.frame() %>% rownames_to_column("genome") %>% data.table() %>% 
  merge(metadata_qc[,list(genome, Continent, Country,  IsMAG, GenomeSource.x)], by = "genome")
pcoa_plot_mash_siso <- ggplot(pcoa4_res, aes(x = Axis.1, y = Axis.2, fill = Continent, shape = ifelse(IsMAG, "MAG", "Isolate"))) + 
  geom_point(size = 1.7) + 
  scale_shape_manual(values = c(21,22), name = "")+ 
  scale_fill_manual(values = continent_colors, name = "Continent", guide = guide_legend(override.aes = list(shape = 21))) + theme_classic2()
ggsave(pcoa_plot_mash_siso, file = "slackia_pcoa_continent.pdf", width = 4.8, height = 3.2)

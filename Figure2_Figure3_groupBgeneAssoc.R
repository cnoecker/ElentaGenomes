## Figure 2, Figure 3 and Figure S4 - gene differences between groups A and B
library(data.table)
library(tidyverse)
library(vegan)
library(ape)
library(ggtree)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
library(ranger)
library(broom)
library(ggrepel)
theme_set(ggpubr::theme_classic2())
metadata_qc <- readRDS("final_qc_metadata.rds")
## QC data for MAGs
el_metadata <- metadata_qc[grepl("lenta", species)]

# Figure 2A-B
genes_plot <- ggplot(el_metadata[!is.na(Clust2) & IsMAG == F], aes(x = factor(Clust2), y = NumCDSBakta)) + 
  geom_violin() + geom_jitter(height = 0, size = 0.5) + theme_classic2() + 
  xlab("Genome Group") + ylab("Number of annotated genes")
#Smaller genomes in Clust2
size_plot <- ggplot(el_metadata[!is.na(Clust2) & IsMAG == F], aes(x = factor(Clust2), y = Genome_Size)) + geom_violin() + geom_jitter(height = 0, size = 0.5) + theme_classic2() + 
  xlab("Genome Group") + ylab("Genome size (bp)")
compare_plots <- size_plot + stat_compare_means() + genes_plot + stat_compare_means()
ggsave(compare_plots, file = "genome_size_genes_clusters_baktaGenes.pdf",
       width = 7, height = 3.5)

# Figure 2C
### Size of core genome in Clust 1 vs Clust 2
elenta_eggnog <- readRDS("elenta_eggnog.rds")
elenta_eggnog <- merge(elenta_eggnog, metadata_qc[,list(genome, Clust2, IsMAG, GenomeSource.x)], by = "genome")
clust_genome_count <- elenta_eggnog[,length(unique(genome)), by = list(Clust2, eggNOG_OGs)]
tot_genomes <- metadata_qc[,length(genome),by=Clust2]
clust_genome_count <- merge(clust_genome_count, tot_genomes, by = "Clust2")
clust_genome_count[,FracGenomes:=V1.x/V1.y]
clust_genome_count[,IsCore:=ifelse(FracGenomes > 0.95, "Core (>95%)", "Accessory")]
clust_pangenome_bar <- ggplot(clust_genome_count, aes(x = factor(Clust2), fill = IsCore)) + geom_bar(color = "black") + 
  scale_fill_brewer(palette = "Purples", name = "") + 
  ylab("Number of eggNOG gene families") + xlab("Genome Group")
ggsave(clust_pangenome_bar, file = "core_pan_cluster.pdf", width = 4.3, height = 3.2)

# Figure 2D
### Differences in gene categories by cluster
top_cat <- elenta_eggnog[,length(unique(eggNOG_OGs)), by=COG_category][order(V1, decreasing = T)][1:10, COG_category]
elenta_eggnog[,COG_category2:=ifelse(COG_category %in% top_cat, COG_category, "Multiple/Other")]
cat_names <- data.table(COG_category = c(elenta_eggnog[nchar(COG_category)==1, sort(unique(COG_category))], "W", "Z"), 
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

clust_genome_count <- elenta_eggnog[,length(unique(genome)), by = list(Clust2, eggNOG_OGs, COG_category, COG_category2)]
tot_genomes <- el_metadata[,length(genome),by=Clust2]
clust_genome_count <- merge(clust_genome_count, tot_genomes, by = "Clust2")
clust_genome_count <- merge(clust_genome_count, cat_names, by.x = "COG_category2", by.y = "COG_category", all.x = T)
clust_genome_count[,FracGenomes:=V1.x/V1.y]
clust_genome_count[COG_category2 %in% c("-", "S", "Multiple/Other"), COG_category2:="Multiple/Other/Unknown"]
clust_genome_count[,COG_category3:=ifelse(COG_category2=="Multiple/Other/Unknown", COG_category2, 
                                          paste0(COG_category2, " ", CategoryName))]

cog_order <- clust_genome_count[FracGenomes > 0.95,length(eggNOG_OGs), by=COG_category3][order(V1, decreasing = T)][,COG_category3]
clust_genome_count[,COG_category3:=factor(COG_category3, levels = c(cog_order))]

#Color scales
cog_colors <- c("gray20", brewer.pal(7, "Dark2"), "gray50", brewer.pal(3, "Set1")[2], 'gray70')
names(cog_colors) <- c(elenta_eggnog[COG_category2 != "Multiple/Other", sort(unique(COG_category2))], "Multiple/Other")

cog_colors2 <- cog_colors[names(cog_colors) %in% elenta_eggnog[,COG_category2]]
cog_colors2 <- c(cog_colors2[c(2:8, 10)], "gray80")
names(cog_colors2) <- c(sort(cog_order[cog_order != "Multiple/Other/Unknown"]), "Multiple/Other/Unknown")
#names(cog_colors2) <- c(clust_genome_count[COG_category3 != "Multiple/Other/Unknown",sort(unique(as.character(COG_category3)))], "Multiple/Other/Unknown")

core_genome_compare <- ggplot(clust_genome_count[FracGenomes > 0.95][order(COG_category3)], aes(x = Clust2, fill = COG_category3)) + 
  geom_bar(position = "fill", color = "black") +
  scale_fill_manual(values = cog_colors2, name = "COG category")+ theme_classic2() + 
  xlab("Genome group") + ylab("Fraction of core genome")
ggsave(core_genome_compare, file = "core_genome_clust_compare.pdf", width = 6.5, height = 5)

# Accessory genome between each? 
accessory_compare <- ggplot(clust_genome_count[FracGenomes < 0.95][order(COG_category3)], aes(x = Clust2, fill = COG_category3)) + 
  geom_bar(position = "fill", color = "black") +
  scale_fill_manual(values = cog_colors2, name = "COG category")+ theme_classic2() + 
  xlab("Genome group") + ylab("Fraction of accessory genome") + scale_x_discrete(labels = c("A", "B"))
clust_genome_count[FracGenomes < 0.95][,fisher.test(COG_category3, Clust2, simulate.p.value = T)]
ggsave(accessory_compare, file = "accessory_genome_clust_compare.pdf", width = 6.5, height = 5)

clust_genome_count[FracGenomes < 0.95, length(unique(eggNOG_OGs)), by = list(Clust2, COG_category3)] %>% dcast(COG_category3~Clust2, value.var = "V1") %>% 
  mutate(Frac1 = `1`/sum(`1`), Frac2 = `2`/sum(`2`))

elenta_wide <- dcast(elenta_eggnog, genome~ eggNOG_OGs, value.var = "Present", fun.aggregate = sum)
prev_genes <- elenta_eggnog[,length(unique(genome)), by = eggNOG_OGs][V1 > 28, eggNOG_OGs] #10%
elenta_wide2 <- elenta_wide[,names(elenta_wide) %in% c("genome", prev_genes, "Clust2"), with = F]
universal <- apply(elenta_wide2[,2:(ncol(elenta_wide2)-1), with = F], 2, function(x){ length(unique(x))==1})
prev_genes <- prev_genes[!universal] #remove genes that are identical
elenta_wide2 <- elenta_wide[,names(elenta_wide) %in% c("genome", prev_genes, "Clust2"), with = F]
prev_genes <- elenta_eggnog[,length(unique(genome)), by = eggNOG_OGs][V1 > 28 & V1 < 265, eggNOG_OGs] # Remove core genes
elenta_wide3 <- elenta_wide[,c("genome", prev_genes), with = F]
elenta_wide3 <- merge(elenta_wide3, el_metadata[,list(genome, Clust2, IsMAG, Genome_Size, Completeness)], by = "genome")
cog_cat_tab <- unique(elenta_eggnog[,list(eggNOG_OGs, TopLevel_eggNOG, COG_category, COG_category2)])

gene_assoc <- rbindlist(lapply(prev_genes, function(x){
  elenta_wide3[,Tmp:=ifelse(get(x)> 0, 1, 0)]
  if(elenta_wide3[,length(unique(Tmp)) > 1]){
    foo <- elenta_wide3[,glm(Tmp~factor(Clust2)+Genome_Size+factor(IsMAG), family = "binomial")] %>% 
      broom::tidy()
    foo <- foo %>% mutate(GeneCluster = x)
    return(foo)
  } else{ return(NULL)}
}))
gene_assoc <- gene_assoc[term != "(Intercept)"]
gene_assoc[,Padj:=p.adjust(p.value, method = "BH"), by = term]
gene_assoc[Padj < 0.1 & grepl("Clust2", term)][order(estimate)]
gene_assoc[grepl("Clust2", term),table(Padj < 0.1)]
gene_assoc[Padj < 0.1 & grepl("Clust2", term)][,table(estimate > 0)] 
gene_assoc[grepl("Clust2", term)][order(p.value)][1:10]

### Volcano plot with associations - Figure 2E
cog_cat_tab <- unique(elenta_eggnog[,list(eggNOG_OGs, TopLevel_eggNOG, COG_category, COG_category2)])
gene_assoc <- merge(gene_assoc, cog_cat_tab, by.x = "GeneCluster", by.y = "eggNOG_OGs")
gene_assoc <- merge(gene_assoc, cat_names, by.x = "COG_category2", by.y = "COG_category", all.x = T)
ggplot(gene_assoc[grepl("Clust2", term)], aes(x = estimate, y = -1*log10(Padj), 
                                              fill = COG_category2)) + geom_point(shape = 21, size = 2) + 
  scale_fill_manual(values = cog_colors) + theme_classic2()

gene_assoc[COG_category2 %in% c("-", "S", "Multiple/Other"), COG_category2:="Multiple/Other/Unknown"]
gene_assoc[,COG_category3:=ifelse(COG_category2=="Multiple/Other/Unknown", COG_category2, 
                                  paste0(COG_category2, " ", CategoryName))]
cog_colors2 <- cog_colors[names(cog_colors) %in% gene_assoc[,COG_category2]]
cog_colors2 <- c(cog_colors2, "gray80")
names(cog_colors2) <- c(gene_assoc[COG_category3 != "Multiple/Other/Unknown",sort(unique(COG_category3))], "Multiple/Other/Unknown")
gene_assoc_sub <- gene_assoc[grepl("Clust2", term)]

cog_assoc_plot <- ggplot(gene_assoc_sub, aes(x = estimate, y = -1*log10(Padj), 
                                             fill = COG_category3)) + geom_hline(yintercept = 1, linetype = 2, color = "gray") + 
  geom_vline(xintercept = 1, linetype = 2, color = "gray") + geom_vline(xintercept = -1, linetype = 2, color = "gray") + 
  geom_point(data = gene_assoc_sub[COG_category2 == "Multiple/Other/Unknown"], shape = 21, size = 2) + 
  geom_point(data = gene_assoc_sub[COG_category2 != "Multiple/Other/Unknown"], shape = 21, size = 2) + 
  scale_fill_manual(values = cog_colors2, name = "COG category") + theme_classic2() + xlim(-7, 7) +
  xlab("Adjusted association with Group 2") 
ggsave(cog_assoc_plot, file = "cog_assoc_clust2.pdf", width = 7.5, height = 3.7)

## Figure 2F - COG1053 and 0243
tmtf_dat <- fread("TFcompare.proteinortho.tsv")
cog1053 <- fread("COG1053compare.proteinortho.tsv")
cog0243 <- fread("COG0243compare.proteinortho.tsv")

setnames(tmtf_dat, gsub(".faa$", "", names(tmtf_dat)))
setnames(cog1053, gsub(".faa$", "", names(cog1053)))
setnames(cog0243, gsub(".faa$", "", names(cog0243)))

tmtf_dat[,GeneCluster:=1:nrow(tmtf_dat)]
tmtf_dat <- melt(tmtf_dat, id.var = c(names(tmtf_dat)[1:3], "GeneCluster"), variable.name = "genome", variable.factor = F)
tmtf_dat[,Present:=ifelse(value == "*", 0, 1)]
tmtf_dat[,CopyNumber:=case_when(
  Present==0 ~ 0, 
  TRUE ~ str_count(value, ",")+1
)]
counts_per_genome <- tmtf_dat[,list(length(unique(GeneCluster[Present==1])), sum(CopyNumber)), by= genome] 

cog1053[,GeneCluster:=1:nrow(cog1053)]
cog1053 <- melt(cog1053, id.var = c(names(cog1053)[1:3], "GeneCluster"), variable.name = "genome", variable.factor = F)
cog1053[,Present:=ifelse(value == "*", 0, 1)]
cog1053[,CopyNumber:=case_when(
  Present==0 ~ 0, 
  TRUE ~ str_count(value, ",")+1
)]
counts_per_genome1053 <- cog1053[,list(length(unique(GeneCluster[Present==1])), sum(CopyNumber)), by= genome] 

cog0243[,GeneCluster:=1:nrow(cog0243)]
cog0243 <- melt(cog0243, id.var = c(names(cog0243)[1:3], "GeneCluster"), variable.name = "genome", variable.factor = F)
cog0243[,Present:=ifelse(value == "*", 0, 1)]
cog0243[,CopyNumber:=case_when(
  Present==0 ~ 0, 
  TRUE ~ str_count(value, ",")+1
)]
counts_per_genome0243 <- cog0243[,list(length(unique(GeneCluster[Present==1])), sum(CopyNumber)), by= genome] 

combined_counts <- merge(counts_per_genome, counts_per_genome1053, by = "genome") %>% 
  merge(counts_per_genome0243)
setnames(combined_counts, c("genome", "TMTFUniq", "TMTFCount", "COG1053Uniq", "COG1053Count", "COG0243Uniq", "COG0243Count"))
combined_counts <- merge(combined_counts, el_metadata[,list(genome, Clust2, Continent, Completeness, IsMAG, GenomeSource.x, Genome_Size)])

prev_clusters <- cog1053[Present==1, length(unique(genome)), by = GeneCluster][V1 > 28 & V1 < 265, GeneCluster] # > 10% prevalence
cog_wide <- dcast(cog1053[GeneCluster %in% prev_clusters], genome~GeneCluster, value.var = "Present")
setnames(cog_wide, c("genome", paste0("col", prev_clusters)))
cog_wide <- cog_wide %>%   merge(el_metadata[,list(genome, Clust2, Completeness, IsMAG, Genome_Size)])

cog_assoc_b <- rbindlist(lapply(prev_clusters, function(x){
  col_id <- paste0("col", x)
  foo <- cog_wide[,glm(get(col_id)~factor(Clust2)+Genome_Size+factor(IsMAG), family = "binomial")] %>% tidy()
  foo <- foo %>% mutate(GeneCluster = gsub("col", "", x))
  return(foo)
}))
cog_assoc_b <- cog_assoc_b[term != "(Intercept)"]
cog_assoc_b[,Padj:=p.adjust(p.value, method = "BH"), by = term]

prev_clusters2 <- cog0243[Present==1, length(unique(genome)), by = GeneCluster][V1 > 28 & V1 < 265, GeneCluster] # > 10% prevalence
cog_wide2 <- dcast(cog0243[GeneCluster %in% prev_clusters2], genome~GeneCluster, value.var = "Present")
setnames(cog_wide2, c("genome", paste0("col", prev_clusters2)))
cog_wide2 <- cog_wide2 %>%   merge(el_metadata[,list(genome, Clust2, Completeness, IsMAG, Genome_Size)])

cog_assoc2b <- rbindlist(lapply(prev_clusters2, function(x){
  col_id <- paste0("col", x)
  foo <- cog_wide2[,glm(get(col_id)~factor(Clust2)+Genome_Size+factor(IsMAG), family = "binomial")] %>% tidy()
  foo <- foo %>% mutate(GeneCluster = gsub("col", "", x))
  return(foo)
}))
cog_assoc2b <- cog_assoc2b[term != "(Intercept)"]
cog_assoc2b[,Padj:=p.adjust(p.value, method = "BH"), by = term]

cog_assoc_combined <- rbind(data.table(cog_assoc_b, Type = "Flavin-dependent"), data.table(cog_assoc2b, Type = "Molydopterin"))

## Assign seed orthologs to clusts
cog0243_long_genes <- cog0243[value != "*",list(GeneCluster, genome, value)]
cog0243_long_genes[,value:=sapply(value, strsplit, split = ",")]
cog0243_long_genes <- unnest(cog0243_long_genes, cols = c(value)) %>% data.table()
cog0243_long_genes <- merge(cog0243_long_genes, elenta_eggnog[,list(`#query`, seed_ortholog, evalue, score)], by.x = "value", by.y = "#query", all.x = T)
#cog0243_long_genes[,length(unique(seed_ortholog)), by = GeneCluster][V1 != 1] 
#cog0243_long_genes[GeneCluster==46]
best_seed_ortholog <- cog0243_long_genes[,list(seed_ortholog[which.max(score)], max(score)), by = GeneCluster]
#best_seed_ortholog[grepl("Elen_0471", V1)]
best_match <- best_seed_ortholog[,GeneCluster[which.max(V2)], by = V1]
setnames(best_match, c("seed_ortholog", "GeneCluster"))

cog1053_long_genes <- cog1053[value != "*",list(GeneCluster, genome, value)]
cog1053_long_genes[,value:=sapply(value, strsplit, split = ",")]
cog1053_long_genes <- unnest(cog1053_long_genes, cols = c(value)) %>% data.table()
cog1053_long_genes <- merge(cog1053_long_genes, elenta_eggnog[,list(`#query`, seed_ortholog, evalue, score)], by.x = "value", by.y = "#query", all.x = T)
#cog1053_long_genes[,length(unique(seed_ortholog)), by = GeneCluster][V1 != 1] 
#cog1053_long_genes[GeneCluster==46]
best_seed_ortholog1053 <- cog1053_long_genes[,list(seed_ortholog[which.max(score)], max(score)), by = GeneCluster]
#best_seed_ortholog1053[grepl("Elen_2529", V1)]
best_match1053 <- best_seed_ortholog1053[,GeneCluster[which.max(V2)], by = V1]
setnames(best_match1053, c("seed_ortholog", "GeneCluster"))

genes_of_interest <- data.table(
  seed_ortholog = c("479437.Elen_2529", "479437.Elen_2528", paste0("479437.Elen_245", 1:3), "479437.Elen_0360", 
                    "479437.Elen_0690", "479437.Elen_1325", "479437.Elen_0471", "479437.Elen_0363", 
                    "479437.Elen_0107", "479437.Elen_0405", "479437.Elen_0507", "479437.Elen_0629", "479437.Elen_0988", 
                    "479437.Elen_2533", "479437.Elen_2653", "479437.Elen_2695", "479437.Elen_3031"), 
  gene_name = c("cgr2", "cgr1", paste0("progestin genes ", 1:3), "3alphaHSD", "3alphaHSD-loxo", "3betaHSD", "dadH", "ber", "urdA", 
                "rrdA", "cadh", "crdA", "crdC", "crdD_cac4", "srdA", "irdA", "hcdh"), 
  substrates = c("cardenolides", "cardenolides", rep("glucocorticoids", 3), "3alpha sterols", "3alpha sterols", "3beta sterols", 
                 "dopamine", "lignans", "urocanate", "resveratrol", "catechin", "cinnamate", "cinnamate", "cinnamate", "shikimate", 
                 "itaconate", "hydrocaffeic acid")
)
best_match_all <- rbind(data.table(best_match, Type = "Molydopterin"), data.table(best_match1053, Type = "Flavin-dependent"))
best_match_all <- merge(best_match_all, genes_of_interest, by = "seed_ortholog", all.x = T)
best_match_all[,GeneCluster:=as.character(GeneCluster)]
cog_assoc_combined <- merge(cog_assoc_combined, best_match_all, by = c("Type", "GeneCluster"), all.x = T)

cog_assoc_combined[,Label:=ifelse(!is.na(gene_name), gene_name, "")]
cog_assoc_combined_clust2 <- cog_assoc_combined[grepl("Clust2", term)]
assoc_plot <- ggplot(cog_assoc_combined_clust2, aes(x = estimate, y = -1*log10(Padj),fill = Type)) + 
  geom_hline(yintercept= 1, linetype = 2, color = "gray") + geom_vline(xintercept = -1, linetype = 2, color = "gray") + geom_vline(xintercept = 1, linetype = 2, color = "gray") + 
  geom_point(size = 2, shape = 21) + theme_classic() + scale_fill_brewer(palette = "Dark2", name = "Enzyme class") + 
  xlab("Adjusted association with Genome Group 2") +
  geom_text_repel(aes(label = Label), color = "black", size = 2.4, max.overlaps = 50, min.segment.length = 0.2) + xlim(-5.5, 5.5)
ggsave(assoc_plot, file = "final_genome_set_2024-07/respiratoryEnzyme_highRes_Clust2Assoc_labels.pdf", width = 5.2, height = 3.7)

cog_assoc_combined_clust2[, length(GeneCluster), by = (estimate > 1 & Padj < 0.1)]
cog_assoc_combined_clust2[, length(GeneCluster), by = (estimate < -1 & Padj < 0.1)]


## Figure 3 A and B
library(gggenes)

#Figure 3B
## HS strain example - HS42-2A
nitro_example <- fread("annots_out/bakta_elentaOnly/HS42-2A.tsv")
surrounding_genes <- nitro_example[grepl("FBGJHN_123[5-8][0-9]", `Locus Tag`)]
surrounding_genes <- nitro_example[Start > 2846000 & Stop < 2854000]
surrounding_genes <- merge(surrounding_genes, elenta_eggnog[,list(`#query`, eggNOG_OGs, Description, PFAMs, TopLevel_eggNOG)], by.x = "Locus Tag", by.y = "#query")
surrounding_genes <- merge(surrounding_genes, gene_assoc_sub, by.x = "eggNOG_OGs", by.y = "GeneCluster", all.x = T)

gene_plot <- ggplot(surrounding_genes, aes(xmin = Stop, xmax = Start, y = `#Sequence Id`, fill = estimate, forward = ifelse(Strand == "-", TRUE, FALSE))) +
  geom_gene_arrow(arrow_body_height = unit(14, "mm"), arrowhead_height = unit(18, "mm")) + #log2(OR+1e-7)
  geom_gene_label(aes(label = str_wrap(Description, width = 14)), color = "white", height = unit(14, "mm")) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = brewer.pal(8, "Blues")[8], midpoint = 0, 
                       name = "Log2 Enrichment in Group B", limits = c(-2, 6)) +
  theme_genes() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
ggsave(gene_plot, file = "nr_region_genes_enrichment2.pdf", width = 10, height = 3.5)

# Figure 3A - Similar gene plot just for Mrr and AbiO region
surrounding_genes <- nitro_example[`Locus Tag` %in% elenta_eggnog[eggNOG_OGs=="COG1715@1|root,COG1715@2|Bacteria,2GP9K@201174|Actinobacteria,4CWR0@84998|Coriobacteriia", `#query`]]
surrounding_genes <- nitro_example[grepl("FBGJHN_105[2-4][0-9]", `Locus Tag`)]
surrounding_genes <- merge(surrounding_genes, elenta_eggnog[,list(`#query`, eggNOG_OGs, Description, PFAMs, TopLevel_eggNOG)], by.x = "Locus Tag", by.y = "#query")
surrounding_genes <- merge(surrounding_genes, gene_assoc_sub, by.x = "eggNOG_OGs", by.y = "GeneCluster", all.x = T)

gene_plot <- ggplot(surrounding_genes, aes(xmin = Stop, xmax = Start, y = `#Sequence Id`, fill = estimate, forward = ifelse(Strand == "-", TRUE, FALSE))) +
  geom_gene_arrow(arrow_body_height = unit(14, "mm"), arrowhead_height = unit(18, "mm")) + #log2(OR+1e-7)
  geom_gene_label(aes(label = str_wrap(Description, width = 14)), color = "white", height = unit(14, "mm")) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = brewer.pal(8, "Blues")[8], midpoint = 0, 
                       name = "log2 Enrichment\nin Group B", limits = c(-2, 6)) +
  theme_genes() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
ggsave(gene_plot, file = "t4re_region_genes_enrichment2.pdf", width = 9, height = 3.5)


## Figure 3C - fast.genomics output
## visualization of fast.genomics results
library(gggenomes)
library(Biostrings)
library(DECIPHER)
library(RColorBrewer)
library(patchwork)

corios <- fread("final_genome_set_2024-07/fastgenomicsOutputNR/browse_260.tsv")
otherbac <- fread("final_genome_set_2024-07/fastgenomicsOutputNR/browse_260(2).tsv")
## Try with gggenomes
setnames(corios, "begin", "start")
corios[,seq_id:=paste0(gtdbSpecies, "_", track)]
top_desc <- corios[,length(gtdbDomain), by = list(clusterId, desc)][,list(desc[which.max(V1)], V1[which.max(V1)]), by = clusterId]
corios <- merge(corios, top_desc, by = "clusterId", all.x = T)
corios[,desc2:=str_trunc(V1, 20)]

nitro_direction <- corios[clusterId=="query", list(track, seq_id, strand, start, end)]
nitro_direction[,seq_id2:="Eggerthella lenta_1"]
nitro_direction[,start2:=161139]
nitro_direction[,end2:=161921]


setnames(otherbac, "begin", "start")
otherbac[,seq_id:=paste0(gtdbSpecies, "_", track)]
top_desc <- otherbac[,length(gtdbDomain), by = list(clusterId, desc)][,list(desc[which.max(V1)], V1[which.max(V1)]), by = clusterId]
otherbac <- merge(otherbac, top_desc, by = "clusterId", all.x = T)
otherbac[,desc2:=str_trunc(V1, 25)]

nitro_direction2 <- otherbac[clusterId=="query", list(track, seq_id, strand, start, end)]


# Need to fix color scales, do tree, add E lenta to otherbac or make a combined version - 
# how to do with tree

seqs1 <- readAAStringSet("fastgenomicsOutputNR/browse_260.fasta")
seqs2 <- readAAStringSet("fastgenomicsOutputNR/browse_260(1).fasta")
seqs_all <- append(seqs1, seqs2)

seqs_all <- seqs_all[-which(duplicated(seqs_all))] # Remove extra Arabiibacter
id_mapping <- data.table(seqName = names(seqs_all))
id_mapping[,seqStrain:=gsub(" assemblyId.*", "", gsub(".*identity from ", "", seqName))]
id_mapping[,seqAssemblyId:=gsub(".* assemblyId ", "", seqName)]
id_mapping[,seqLocusTag:=gsub(" .*", "", seqName)]
id_mapping <- merge(id_mapping, unique(corios[,list(assemblyId, locusTag, seq_id)]), 
                    by.x = c("seqAssemblyId", "seqLocusTag"), by.y = c("assemblyId", 
                                                                       "locusTag"), all.x = T)
id_mapping <- merge(id_mapping, unique(otherbac[,list(assemblyId, locusTag, seq_id)]), 
                    by.x = c("seqAssemblyId", "seqLocusTag"), by.y = c("assemblyId", 
                                                                       "locusTag"), all.x = T)
id_mapping[,seq_id:=ifelse(is.na(seq_id.x), seq_id.y, seq_id.x)]
id_mapping[,CoriosOther:=ifelse(is.na(seq_id.x), "Other", "Corios")]
id_mapping[,bitScore:=as.numeric(gsub(".* ", "", gsub(" bits.*", "", seqName)))]
id_mapping[,percIdent:=as.numeric(gsub(".* ", "", gsub("% identity.*", "", seqName)))]
#id_mapping[order(percIdent, decreasing = T)]

seqs_all_sub <- seqs_all[names(seqs_all) %in% id_mapping[order(percIdent, decreasing = T)][1:25, seqName]]
seqs_all_sub <- AlignSeqs(seqs_all_sub, gapOpening = -55, gapExtension = -40, iterations = 3, refinements = 5)

names(seqs_all_sub) <- id_mapping[match(names(seqs_all_sub), seqName), seq_id]

seqs_tree <- Treeline(seqs_all_sub)

combined_gene_dat <- rbind(corios, otherbac)
combined_gene_dat <- combined_gene_dat[seq_id %in% names(seqs_all_sub)]
combined_gene_dat[clusterId=="query" & grepl("lenta", gtdbSpecies), strand]
## Next, figure out color scale

combined_gene_dat[,desc3:=case_when(
  grepl("nitroreductase", desc) ~ "nitroreductase family", 
  grepl("EFR1", desc) ~ "EFR1 family ferrodoxin", 
  desc == "TetR/AcrR family transcriptional regulator" ~ "TetR family regulator",
  desc == "glutathione peroxidase" ~ "glutathione peroxidase",
  TRUE ~ clusterId
)]
combined_gene_dat[,desc3:=factor(desc3, levels = c(
  "", c(1:42)[1:42 %in% combined_gene_dat[,desc3]], 
  "nitroreductase family",
  "EFR1 family ferrodoxin",
  "TetR family regulator",
  "glutathione peroxidase"
))]
color_scale1 <- c(paste0("gray", seq(10, 83, by = 3)), brewer.pal(4, "Set1"))
names(color_scale1) <- levels(combined_gene_dat[,desc3])

combined_gene_dat_sub <- combined_gene_dat %>% filter(
  !grepl("urolithinfacients_7", seq_id) & 
    !grepl("urolithinfacients_2", seq_id) & 
    !grepl("rubneri_1", seq_id) & 
    !grepl("Rubneribacter", seq_id) & 
    !grepl("lentibacter hominis", seq_id) & 
    !grepl("Beduin", seq_id) & 
    !grepl("Eubacterium", seq_id) &
    !grepl("sp000270285", seq_id)
)

nitro_sub <- combined_gene_dat_sub[grepl("nitroreductase", desc3)]
combined_gene_dat_sub <- merge(combined_gene_dat_sub, 
                               nitro_sub[,list(seq_id, start, end)], by = "seq_id", all.x = T)

combined_gene_dat_sub[,diffStart:=abs(start.x-start.y)]
combined_gene_dat_sub[,diffEnd:=abs(start.x-end.y)]
combined_gene_dat_sub[,diffStart2:=abs(end.x-start.y)]
combined_gene_dat_sub[,diffEnd2:=abs(end.x-end.y)]
combined_gene_dat_sub[,minDiff:=ifelse(diffStart < diffEnd, diffStart, diffEnd)]
combined_gene_dat_sub[,minDiffEnd:=ifelse(diffStart2 < diffEnd2, diffStart2, diffEnd2)]
combined_gene_dat_sub[,minDiffAll:=ifelse(minDiff < minDiffEnd, minDiff, minDiffEnd)]

combined_gene_dat_sub <- combined_gene_dat_sub %>% filter(
  !(minDiffAll > 4000 & (desc3 == ""|is.na(desc3)))
)
combined_gene_dat_sub[,start:=start.x]
combined_gene_dat_sub[,end:=end.x]

seqs_all_sub2 <- seqs_all_sub[names(seqs_all_sub) %in% combined_gene_dat_sub[,seq_id]]
seqs_tree2 <- Treeline(seqs_all_sub2)

tree_plot2 <- ggtree::ggtree(seqs_tree2, branch.length = "none") + 
  ggtree::geom_tiplab() + ggtree::geom_nodelab()

tree_plot_order2 <- tree_plot2$data %>% select(label, y) %>% filter(!is.na(label)) %>% 
  arrange(y) %>% pull(label)

combined_gene_dat_sub[,seq_id:=factor(seq_id, labels = rev(tree_plot_order2))]
gene_plot2 <- gggenomes(genes = combined_gene_dat_sub) %>% 
  flip(combined_gene_dat_sub[clusterId=="query" & strand == "+", seq_id]) + 
  geom_seq() + 
  geom_gene(aes(fill = desc3)) +
  geom_seq_label(aes(label = seq_id), size = 4) +
  scale_fill_manual(values = color_scale1, 
                    breaks = names(color_scale1)[26:29], name = "", 
                    labels = c("nitroreductase family", 
                               "4Fe-4S protein", 
                               "TetR family regulator", 
                               "glutathione peroxidase"))

combined_plot <- tree_plot2+gene_plot2+ plot_layout(widths = c(1, 4))
ggsave(combined_plot, file = "nitroreductase_fast_genomics_viz.pdf", width = 10, height = 7.5)

## Figure S4 - gene regions
# Save top hits

# Random forest to predict Group B
rf1 <- ranger(dependent.variable.name = "Clust2", data = elenta_wide2[,2:ncol(elenta_wide2), with=F],
              importance = "permutation", probability = T)

probabilities <- data.frame(pred = predict(rf1, data = elenta_wide2[,2:(ncol(elenta_wide2)-1), with=F], 
                                           type='response', verbose = TRUE)$predictions)
probabilities$binary = ifelse(probabilities[,2] > probabilities[,1], 2, 1)
roc1 <- roc(response = el_metadata[, Clust2], predictor = probabilities[,2])
table(probabilities$binary, elenta_wide2[,Clust2])
auc(roc1)
# Save importance scores for each gene family
imp_results <- data.table(GeneFam = names(rf1$variable.importance), Importance = rf1$variable.importance)

# gene_assoc is from Fig 2 above
genes_combined_tab <- merge(gene_assoc[grepl("Clust2", term)], imp_results, by.x = "GeneCluster", by.y = "GeneFam", all = T)
top_genes_combined_tab <- genes_combined_tab[estimate > 0 & ((estimate > 2 & Padj < 1e-10) | Importance > 0.005)]  
top_genes_combined_tab <- top_genes_combined_tab[order(estimate), list(GeneCluster, estimate, p.value, Padj, Importance, COG_category2, Description, Taxon)]

genes_sub <- elenta_eggnog[eggNOG_OGs %in% top_genes_combined_tab[,GeneCluster]]
genomes_sub <- genes_sub[,sort(unique(genome))]

# Is there a genome with all the top genes?
elenta_eggnog[eggNOG_OGs %in% top_genes_combined_tab[,GeneCluster], length(unique(eggNOG_OGs)), by = genome][order(V1, decreasing = T)][1:10]

top_imp_genes <- top_genes_combined_tab[Importance > 0.004 | estimate > 4.5]
genome_example <- fread("annots_out/bakta_elentaOnly/GCF_021378625.1_ASM2137862v1_genomic.tsv")

genome_example_hits <- genome_example[`Locus Tag` %in% elenta_eggnog[eggNOG_OGs %in% top_imp_genes[,GeneCluster], `#query`]]
regions <- data.table(
  start = c(150000, 2225000, 2326000, 2438000, 2748000), end = c(205000, 2232000, 2365000, 2464000, 2766000)
)
genome_example_regions <- genome_example[
  (Start > regions[1,start] & Stop < regions[1,end]) |
    (Start > regions[2,start] & Stop < regions[2,end]) |
    (Start > regions[3,start] & Stop < regions[3,end]) |
    (Start > regions[4,start] & Stop < regions[4,end]) |
    (Start > regions[5,start] & Stop < regions[5,end])
]
genome_example_regions <- merge(genome_example_regions, elenta_eggnog[,list(`#query`, eggNOG_OGs, Description, PFAMs, TopLevel_eggNOG)], by.x = "Locus Tag", by.y = "#query")
genome_example_regions <- merge(genome_example_regions, genes_combined_tab[,list(GeneCluster, estimate, Padj, Importance)], 
                                by.x = "eggNOG_OGs", by.y = "GeneCluster", all.x = T)
genome_example_regions <- genome_example_regions[order(Start)]

gene_plot1 <- ggplot(genome_example_regions[1:36], aes(xmin = Stop, xmax = Start, y = `#Sequence Id`, fill = estimate, 
                                                       forward = ifelse(Strand == "-", TRUE, FALSE))) +
  geom_gene_arrow(arrow_body_height = unit(14, "mm"), arrowhead_height = unit(18, "mm")) + 
  geom_gene_label(aes(label = str_wrap(Description, width = 14)), color = "white", height = unit(14, "mm")) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = brewer.pal(8, "Blues")[8], midpoint = 0, limits = c(-2, 7),
                       name = "log2 Enrichment in Group 2") +
  theme_genes() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
ggsave(gene_plot1, file = "group2_region1_enrichment.pdf", width = 25, height = 3.5)

gene_plot2 <- ggplot(genome_example_regions[Start > regions[2, start] & Stop < regions[2, end]], aes(xmin = Stop, xmax = Start, 
                                                                                                     y = `#Sequence Id`, 
                                                                                                     fill = estimate, 
                                                                                                     forward = ifelse(Strand == "-", TRUE, FALSE))) +
  geom_gene_arrow(arrow_body_height = unit(14, "mm"), arrowhead_height = unit(18, "mm")) + 
  geom_gene_label(aes(label = str_wrap(Description, width = 14)), color = "white", height = unit(14, "mm")) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = brewer.pal(8, "Blues")[8], midpoint = 0, limits = c(-2, 7),
                       name = "log2 Enrichment in Group 2") +
  theme_genes() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
ggsave(gene_plot2, file = "group2_region2_enrichment.pdf", width = 8, height = 3.5)

gene_plot3 <- ggplot(genome_example_regions[Start > regions[3, start] & Stop < regions[3, end]], aes(xmin = Stop, xmax = Start, 
                                                                                                     y = `#Sequence Id`, 
                                                                                                     fill = estimate, 
                                                                                                     forward = ifelse(Strand == "-", TRUE, FALSE))) +
  geom_gene_arrow(arrow_body_height = unit(14, "mm"), arrowhead_height = unit(18, "mm")) + 
  geom_gene_label(aes(label = str_wrap(Description, width = 14)), color = "white", height = unit(14, "mm")) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = brewer.pal(8, "Blues")[8], midpoint = 0, limits = c(-2, 7),
                       name = "log2 Enrichment in Group 2") +
  theme_genes() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
ggsave(gene_plot3, file = "group2_region3_enrichment.pdf", width = 24, height = 3.5)

gene_plot4 <- ggplot(genome_example_regions[Start > regions[4, start] & Stop < regions[4, end]], aes(xmin = Stop, xmax = Start, 
                                                                                                     y = `#Sequence Id`, 
                                                                                                     fill = estimate, 
                                                                                                     forward = ifelse(Strand == "-", TRUE, FALSE))) +
  geom_gene_arrow(arrow_body_height = unit(14, "mm"), arrowhead_height = unit(18, "mm")) + 
  geom_gene_label(aes(label = str_wrap(Description, width = 14)), color = "white", height = unit(14, "mm")) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = brewer.pal(8, "Blues")[8], midpoint = 0, limits = c(-2, 7),
                       name = "log2 Enrichment in Group 2") +
  theme_genes() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
ggsave(gene_plot4, file = "group2_region4_enrichment.pdf", width = 12, height = 3.5)

gene_plot5 <- ggplot(genome_example_regions[Start > regions[5, start] & Stop < regions[5, end]], aes(xmin = Stop, xmax = Start, 
                                                                                                     y = `#Sequence Id`, 
                                                                                                     fill = estimate, 
                                                                                                     forward = ifelse(Strand == "-", TRUE, FALSE))) +
  geom_gene_arrow(arrow_body_height = unit(14, "mm"), arrowhead_height = unit(18, "mm")) + 
  geom_gene_label(aes(label = str_wrap(Description, width = 14)), color = "white", height = unit(14, "mm")) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = brewer.pal(8, "Blues")[8], midpoint = 0, limits = c(-2, 7),
                       name = "log2 Enrichment in Group 2") +
  theme_genes() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
ggsave(gene_plot5, file = "group2_region5_enrichment.pdf", width = 19, height = 3.5)


## Save tables with top genes and annots
top_genes_annots_sub <- merge(elenta_eggnog[eggNOG_OGs %in% top_genes_combined_tab[,GeneCluster]], 
                              bakta_dat_genes_sub[,list(`Locus Tag`, Gene, Product, DbXrefs)], 
                              by.x = "#query", by.y = "Locus Tag", all.x = T, all.y = F)
top_genes_annots_sub <- top_genes_annots_sub[,list(eggNOG_OGs, Description, 
                                                   PFAMs, Gene, Product, DbXrefs)] %>% distinct()
top_genes_annots_sub[,Uniref90:=gsub(".*UniRef90_", "", DbXrefs)]
fwrite(top_genes_annots_sub, file = "top_genes_annots_eggNOG_Uniref.tsv", sep = "\t")

# Get Group1-assoc genes also, for humann comparison
group1_geneclusts <- merge(
  elenta_eggnog[eggNOG_OGs %in% group1_genes[,GeneCluster]], bakta_dat_genes_sub[,list(`Locus Tag`, Product, DbXrefs)],  
  by.x = "#query", by.y = "Locus Tag", all.x = T, all.y = F
)
group1_geneclusts <- group1_geneclusts[,list(eggNOG_OGs, Description, Product, DbXrefs)] %>% distinct()
group1_geneclusts[,Uniref90:=gsub(".*UniRef90_", "", DbXrefs)]
group1_geneclusts[Uniref90==""|grepl(",", Uniref90)|grepl("PFAM:", Uniref90), Uniref90:=NA]
group1_geneclusts <- group1_geneclusts[!is.na(Uniref90)]
group1_geneclusts[,length(unique(Uniref90)), by = eggNOG_OGs][,hist(V1, breaks = 10)]
group1_geneclusts <- merge(group1_geneclusts, group1_genes[,list(GeneCluster, estimate, statistic, p.value, Padj)], 
                           by.x = "eggNOG_OGs", by.y = "GeneCluster", all.x = T)
fwrite(group1_geneclusts, file = "group1_assoc_genes_eggNOG.tsv", sep = "\t")

core_genes_tab <- elenta_eggnog[SoftCore==1]
core_genes_tab <- merge(
  elenta_eggnog[SoftCore==1], bakta_dat_genes_sub[,list(`Locus Tag`, Product, DbXrefs)],  
  by.x = "#query", by.y = "Locus Tag", all.x = T, all.y = F
)
core_genes_tab <- core_genes_tab[,list(eggNOG_OGs, Description, Product, DbXrefs)] %>% distinct()
core_genes_tab[,Uniref90:=gsub(".*UniRef90_", "", DbXrefs)]
core_genes_tab[Uniref90==""|grepl(",", Uniref90)|grepl("PFAM:", Uniref90), Uniref90:=NA]
core_genes_tab <- core_genes_tab[!is.na(Uniref90)]
core_genes_tab[,length(unique(Uniref90)), by = eggNOG_OGs][,hist(V1, breaks = 50)]
fwrite(core_genes_tab, file = "core_genes_eggnog_Uniref.tsv", sep = "\t")

## Figure 3D - NO tolerance data
library(data.table)
library(readxl)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

theme_set(theme_classic2())
# Import growth curve analysis functions
source("growth_curve_functions_byonoy.R")

datadir <- "nitro_reanalysis//"
## Data files could be all csv files in the datadir
datafiles <- list.files(datadir, pattern = ".csv", full.names = T)

metadata_files <- sort(datafiles[grepl("Layout", datafiles)])
datafiles <- sort(datafiles[!datafiles %in% metadata_files])

# Read in all the metadata files
met1 <-  read_csv(metadata_files[1])
met1 <- mutate(met1, Experiment = gsub("PlateMap.*", "", basename(metadata_files[1])))
met1 <- data.table(met1)
met1 <- met1[1:16]
met1[,Type:=rep(c("Condition", "Strain"), 8)]
setnames(met1, "...1", "Row")
met1 <- fread(metadata_files[1], skip = 7, header = F)
met1 <- rbind(met1, met1)
met1 <- mutate(met1, Type = c(rep("Condition", 8), rep("Strain", 8)))
met1[,V2:=gsub("No NO", "NoNO", V2)]
met1[,V13:=gsub("No NO", "NoNO", V13)]
met1[,V10:=gsub("UCSF 2243", "UCSF2243", V10)]
met1[,V11:=gsub("UCSF 2243", "UCSF2243", V11)]
for(j in 2:13){
  foo <- names(met1)[j]
  met1[,c(j):=ifelse(Type != "Condition", gsub(".* ", "", get(..foo)),
                     gsub(" .*", "", get(..foo)))]
}
setnames(met1, names(met1)[1:13], c("Row", 1:12)) 
met1[,Experiment:=gsub("PlateMap.*", "", basename(metadata_files[1]))]
met2 <- fread(metadata_files[2],skip = 1, header = F)
met2 <- met2[1:16]
met2[,Experiment:=gsub("Layout.csv", "", basename(metadata_files[2]))]
met2[,Type:=rep(c("Condition", "Strain"), 8)]
#met3 <- mutate(met3, Type = c(rep("Condition", 8), rep("Strain", 8)))
setnames(met2, names(met2)[1:13], names(met1)[1:13])

met3 <- fread(metadata_files[3],skip = 1, header = F)
met3 <- met3[1:16]
met3[,Experiment:=gsub("Layout.csv", "", basename(metadata_files[3]))]
met3[,Type:=rep(c("Condition", "Strain"), 8)]
setnames(met3, names(met3)[1:13], names(met1)[1:13])
met4 <- fread(metadata_files[4], skip = 1, header = F)
met4[,Experiment:=gsub("Layout.*", "", basename(metadata_files[4]))]
met4 <- met4[1:16]
met4[,Type:=rep(c("Condition", "Strain"), 8)]
setnames(met4, names(met4)[1:13], names(met1)[1:13])

## Convert to long format
layout_info <- rbind(met1, met2, met3, met4)
layout_info <- layout_info %>% pivot_longer(-any_of(c("Row", "Type", "Experiment")))
layout_info <- pivot_wider(layout_info, names_from = "Type", values_from = "value")
names(layout_info)[3] <- "Column"
layout_info <- layout_info %>% 
  mutate(well = paste0(Row, Column))
layout_info <- layout_info %>% 
  mutate(Condition = gsub("um", "uM", (gsub("No NO", "NoNO", Condition))))
## Make condition into an ordered set of categories
layout_info$Condition <- factor(layout_info$Condition, levels = c("NoNO", "16uM", "32uM", "64uM","128uM", "256uM"))
layout_info <- layout_info %>% 
  mutate(Strain=gsub("-[1-3]$", "", Strain))
layout_info <- layout_info %>% 
  mutate(Strain = ifelse(grepl("2243", Strain), "DSM2243", Strain))
layout_info <- layout_info %>% 
  mutate(Strain = ifelse(grepl("HS111", Strain), "HS111-1A", Strain))
#For normalization, the metadata needs to have a column named "Condition" to specify which blank wells/controls 

# Get set of blank wells for normalization from metadata
blank_wells1 = layout_info %>% filter(Strain=="CTRL") %>% pull(well)
names(blank_wells1) = layout_info %>% filter(Strain == "CTRL") %>% pull(Condition)

# Read in all the metadata files
growth_data1 <- data.frame() # An empty table that will hold all the metadata
for(i in 1:length(datafiles)){ #Go through each growth curve file and do this thing
  gtab <- read_growthcurve_file_byonoy(datafiles[i])
  gtab <- mutate(gtab, Experiment = basename(datafiles[i]))
  growth_data1 <- rbind(growth_data1, gtab)
}
growth_data1 <- mutate(growth_data1, Experiment = 
                         case_when(
                           grepl("PILOT3Redo", Experiment) ~ "Pilot 3 Redo",
                           grepl("[PILOT3] Results", Experiment, fixed = T) ~ "Pilot 3",
                           grepl("Pilot6", Experiment) ~ "Pilot 6",
                           grepl("Pilot7", Experiment) ~ "Pilot 7"
                         ))
layout_info <- mutate(layout_info, Experiment = 
                        case_when(
                          grepl("PILOT3Redo", Experiment) ~ "Pilot 3 Redo",
                          grepl("Pilot3", Experiment) ~ "Pilot 3",
                          grepl("Pilot6", Experiment) ~ "Pilot 6",
                          grepl("Pilot7", Experiment) ~ "Pilot 7"
                        ))
# Formatting and initial plot
growth_data_long = pivot_longer(growth_data1, cols = A1:H12, names_to = "well")
growth_data_long = left_join(growth_data_long, layout_info, by = join_by(Experiment, well))
growth_data_long = growth_data_long %>% filter(!is.na(Condition)) #wells without metadata were empty

# Normalize based on early time points
experiments <- c("Pilot 3", "Pilot 3 Redo", "Pilot 6", "Pilot 7")
corrected_data <- data.frame()
for(j in experiments){
  foo <- filter(growth_data_long, Experiment == j)
  corrected_data1 = correct_growth_data_custom(foo, 
                                               method = "earlyData", ctrlMinTime = 0.5, 
                                               ctrlMaxTime = 2, wide_form = F)
  corrected_data <- rbind(corrected_data, corrected_data1)
}
#Set up normalized data
growth_data_norm = dcast(corrected_data, Experiment+time~well, value.var = "NewValue", fill = 0)

# Exclude abnormal readings
layout_info <- layout_info %>% filter(!(Experiment=="Pilot 7" & Strain == "HS23-3A"& well %in% c("D4", "H4")))

el_metadata_sub <- el_metadata[,list(genome, Clust2)]
el_metadata_sub <- el_metadata_sub[genome %in% layout_info[,Strain]|grepl("00024265", genome)|genome=="HS111-1A_CN11"]
el_metadata_sub[,genome:=case_when(
  grepl("HS111", genome) ~ "HS111-1A",
  grepl("24265", genome) ~ "DSM2243",
  TRUE~genome
)]

growth_data_norm = as.data.table(growth_data_norm)
layout_info <-  as.data.table(layout_info)
layout_info <- merge(layout_info, el_metadata_sub, by.x = "Strain", by.y = "genome", all.x = T)
growth_data_norm1 = melt(growth_data_norm, id.var = c("time", "Experiment"), variable.name = "well", variable.factor = F)
growth_data_norm1 = merge(growth_data_norm1, layout_info, by = c("Experiment", "well"), all.x = T)
plot_data = growth_data_norm1
plot1 = ggplot(plot_data[Strain != "" & !is.na(Strain) & 
                           Strain != "CTRL" & Experiment %in% c("Pilot 3", "Pilot 7") & Condition %in% c("NoNO", "64uM", "128uM")], 
               aes(x = time, y = value)) + 
  geom_line(aes(color = Condition, group = well)) + 
  scale_color_manual(values = c("gray30", brewer.pal(3, "Reds")[2:3])) +
  facet_wrap(~paste0(ifelse(Clust2==1, "A", "B"),"_",Strain), scales = "free_y", ncol = 5) +
  ylab("Normalized OD600") + xlab("Time (hr)")

ggsave(plot1, file = "cyrille_NO_growthCurves.pdf", width = 11, height = 3.7)

# Get relative max ODs
corrected_data <- corrected_data %>% filter(!(Experiment=="Pilot 7" & Strain == "HS23-3A"& well %in% c("D4", "H4")))
max_ods <- corrected_data %>% group_by(Experiment, Strain, Condition, well) %>% 
  # Only look for the max value in the latter part of hte growth curve
  filter(time > 0.8*max(time)) %>% 
  summarize(maxOD = max(NewValue))

ctrl_maxes <- max_ods %>% filter(Condition == "NoNO") %>% 
  group_by(Experiment, Strain) %>% summarize(max_ctrl = median(maxOD))

max_ods_compare <- filter(max_ods, Condition != "NoNO") %>% 
  left_join(ctrl_maxes, by = join_by(Strain, Experiment))
max_ods_frac <- max_ods_compare %>% group_by(Experiment, Strain, Condition) %>% 
  summarize(od_frac = median(maxOD, na.rm = T)/mean(max_ctrl))

max_ods_compare %>% 
  filter(Condition %in% c("64uM", "128uM") & Experiment %in% c("Pilot 3", "Pilot 7") & Strain != "CTRL") %>% 
  ggplot(aes(x = Condition, y = maxOD/max_ctrl, fill = Condition)) + 
  stat_summary(geom = "bar", fun = "median") +
  geom_point(shape = 21) + 
  facet_wrap(~Strain) + scale_fill_brewer(palette = "Blues")+ 
  stat_compare_means(method = "anova")+ ylab("Fraction of max OD in no-NO ctrl") +
  ylim(0, 1.2)

max_ods_compare <- left_join(max_ods_compare, el_metadata_sub, by = join_by(Strain == genome))

od_compare_plots <- max_ods_compare %>% 
  filter(Condition %in% c("64uM", "128uM") & Experiment %in% c("Pilot 3", "Pilot 7") & Strain != "CTRL") %>% 
  group_by(Condition, Strain, Clust2) %>%
  summarize(medianODfrac = median(maxOD/max_ctrl)) %>% 
  ggplot(aes(x = ifelse(Clust2==1, "A", "B"), y = medianODfrac)) +
  geom_hline(yintercept = 1, linetype = 2, color = "gray") +
  geom_boxplot(aes(color = factor(Clust2)), outlier.shape = NA)+ 
  geom_jitter(height = 0, width = 0.1) + facet_wrap(~Condition) + 
  scale_color_manual(values = brewer.pal(5,"Blues")[c(3,5)], guide = "none") + theme_classic2() + stat_compare_means() + xlab("Genome Group") +
  ylab("Median max absorbance relative to control") + ylim(-0.02, 1.2) + theme(strip.text = element_text(size = 12))
ggsave(od_compare_plots, file = "maxODcompare_pilots_CN.pdf", width = 4.5, height = 4) 




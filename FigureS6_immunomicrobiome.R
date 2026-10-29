## Figures S6 and S7 - metagenome analysis of Group A/Group B
library(data.table) 
library(tidyverse)
library(Biostrings)
library(readxl)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(ggrepel)
theme_set(theme_classic2())

## Dataset 1: Immunomicrobiome Metagenomic Data - Figure S6
metadata_qc <- readRDS("final_qc_metadata.rds")
el_metadata <- metadata_qc[grepl("lenta", species)]
priv_genes_all <- fread("snps_clust_markers_core_genes.tsv")

source("process_pileup_func.R")
samp_list <- list.files("im_align_out/", pattern = "pileup_out.txt", full.names = T)

pileup_results <- rbindlist(lapply(samp_list, process_pileup, cov_threshold = 0)) # Was 4 previously and only counting major allele count
pileup_results[,Clust2Major:=ifelse(ClusterMarker == "Clust2" & MajorAllele == SNPMarker & MajorAlleleCount > 4, 1, 0)]
pileup_results[,Clust2Minor:=ifelse(ClusterMarker == "Clust2" & MinorAllele == SNPMarker& MinorAlleleCount > 4, 1, 0)]
pileup_results[,Clust1Major:=ifelse(ClusterMarker == "Clust1" & MajorAllele == SNPMarker & MajorAlleleCount > 4, 1, 0)]
pileup_results[,Clust1Minor:=ifelse(ClusterMarker == "Clust1" & MinorAllele == SNPMarker & MinorAlleleCount > 4, 1, 0)]

pileup_results[,Clust2Coverage:=case_when(
  Clust2Major == 1 ~ MajorAlleleCount, 
  Clust2Minor == 1 ~ MinorAlleleCount,
  TRUE ~ 0
)]
pileup_results[,Clust1Coverage:=case_when(
  Clust1Major == 1 ~ MajorAlleleCount, 
  Clust1Minor == 1 ~ MinorAlleleCount,
  TRUE ~ 0
)]
pileup_results_summarize <- pileup_results[,list(
  Clust2MarkersMajor = length(MajorAlleleCount[ClusterMarker == "Clust2" & MajorAlleleCount > 4 & MajorAllele == SNPMarker]),
  Clust2MarkersMinor = length(MinorAlleleCount[ClusterMarker == "Clust2" & MinorAlleleCount > 4 & MinorAllele == SNPMarker]),
  Clust2MarkersMajor2 = length(MajorAlleleCount[ClusterMarker == "Clust2" & MajorAlleleCount > 3 & MajorAllele == SNPMarker]),
  Clust2MarkersMinor2 = length(MinorAlleleCount[ClusterMarker == "Clust2" & MinorAlleleCount > 3 & MinorAllele == SNPMarker]),
  Clust2MedianCount = median(Clust2Coverage[Clust2Major==1 |Clust2Minor==1], na.rm = T),
  Clust2Q25Count = quantile(Clust2Coverage[Clust2Major==1 |Clust2Minor==1], na.rm = T, probs = 0.25),
  Clust2Q75Count = quantile(Clust2Coverage[Clust2Major==1 |Clust2Minor==1], na.rm = T, probs = 0.75),
  Clust2TM = mean(Clust2Coverage[Clust2Major == 1| Clust2Minor == 1], na.rm = T, trim = 0.25),
  Clust2Coverage = sum(Clust2Coverage[Clust2Major==1 |Clust2Minor==1])/length(Clust2Coverage[Clust2Major==1 |Clust2Minor==1]),
  Clust1MarkersMajor = length(MajorAlleleCount[ClusterMarker == "Clust1" & MajorAlleleCount > 4 & MajorAllele == SNPMarker]),
  Clust1MarkersMinor = length(MinorAlleleCount[ClusterMarker == "Clust1" & MinorAlleleCount > 4 & MinorAllele == SNPMarker]),
  Clust1MarkersMajor2 = length(MajorAlleleCount[ClusterMarker == "Clust1" & MajorAlleleCount > 3 & MajorAllele == SNPMarker]),
  Clust1MarkersMinor2 = length(MinorAlleleCount[ClusterMarker == "Clust1" & MinorAlleleCount > 3 & MinorAllele == SNPMarker]),
  Clust1MedianCount = median(Clust1Coverage[Clust1Major==1 | Clust1Minor == 1], na.rm = T),
  Clust1Q25Count = quantile(Clust1Coverage[Clust1Major==1 |Clust1Minor==1], na.rm = T, probs = 0.25),
  Clust1Q75Count = quantile(Clust1Coverage[Clust1Major==1 |Clust1Minor==1], na.rm = T, probs = 0.75),
  Clust1TM = mean(Clust1Coverage[Clust1Major == 1| Clust1Minor == 1], na.rm = T, trim = 0.25),
  Clust1Coverage = sum(Clust1Coverage[Clust1Major==1 | Clust1Minor == 1])/length(Clust1Coverage[Clust1Major==1 | Clust1Minor == 1]),
  TotalCovMedian = as.double(median(V4, trim = 0.25, na.rm = T)), 
  TotalCovMean = as.double(mean(V4, trim = 0.25, na.rm = T)),
  TotMarkersWCoverage = length(Gene[NumCharSimplified > 4]),
  TotMarkersWCoverage2 = length(Gene[NumCharSimplified > 3]),
  TotGenesWCoverage = length(unique(Gene[NumCharSimplified > 4])), 
  TotGenesWCoverage2 = length(unique(Gene[NumCharSimplified > 3]))
), by = Sample]

pileup_results_summarize[,Clust1IQR:=Clust1Q75Count-Clust1Q25Count]
pileup_results_summarize[,Clust2IQR:=Clust2Q75Count-Clust2Q25Count]
readcounts <- fread("readstats.txt", header = F)
pileup_results_summarize <- merge(pileup_results_summarize, readcounts, by.x = "Sample", by.y = "V1", all.x = T)

# IQR cutoff should be % of mean
priv_genes_all[ClustPrevalence > 0.1, length(unique(Gene))] #117 genes, 1224 markers
priv_genes_all[ClustPrevalence > 0.1, length(unique(Gene)), by = ClusterMarker] #75 and 44
priv_genes_all[ClustPrevalence > 0.1, length(Gene)]
pileup_results[,length(unique(Gene))] # only 57 genes
pileup_results[,length(unique(paste0(Gene, Position))), by = ClusterMarker] # only 630 markers

pileup_results_summarize[,Clust1Present:=ifelse(Clust1MarkersMajor + Clust1MarkersMinor >= 8, 1, 0)]
pileup_results_summarize[,Clust2Present:=ifelse(Clust2MarkersMajor + Clust2MarkersMinor >= 8, 1, 0)]

pileup_results_summarize[,Clust1Present2:=ifelse(Clust1MarkersMajor2 + Clust1MarkersMinor2 >= 8, 1, 0)]
pileup_results_summarize[,Clust2Present2:=ifelse(Clust2MarkersMajor2 + Clust2MarkersMinor2 >= 8, 1, 0)]

# if there is not sufficient coverage for at least 20% of genes/markers
pileup_results_summarize[,InsufficientEvidence:=ifelse(TotGenesWCoverage < 21, 1, 0)]

## Import Metaphlan data to compare
metaphlan_dat <- rbindlist(lapply(list.files("metaphlan/"), function(x){
  foo <- fread(paste0("metaphlan/", x))
  el_abund <- foo[grepl("s__Eggerthella_lenta$", `#clade_name`)]
  if(nrow(el_abund) > 0){
    return(data.table(Sample = gsub(".metaphlan", "", x), MetaphlanElenta = el_abund[,relative_abundance], MetaphlanElentaCov = el_abund[,coverage]))
  } else {
    return(data.table(Sample = gsub(".metaphlan", "", x), MetaphlanElenta = 0, MetaphlanElentaCov = 0))
  }
}))
isolate_dat <- metadata_qc[!is.na(ID2),list(ID2, genome, Clust2)]
isolate_dat2 <- fread("extra_hs_genomes_clust_assignment.tsv")
isolate_dat <- rbind(isolate_dat, isolate_dat2[,list(SampleID, ExtraGenome, Clust2)], use.names = F)
isolates_summary <- isolate_dat[,list(NumClust1 = sum(Clust2==1, na.rm = T), NumClust2 = sum(Clust2 == 2, na.rm = T)), by = ID2]
isolates_summary <- merge(isolates_summary, metaphlan_dat, by.x = "ID2", by.y = "Sample", all.x = T)

# Number of markers per sample
clust2_markers_mgs <- ggplot(pileup_results_summarize[!is.na(NumClust2)], aes(x = NumClust2 > 0,
                                                                              y = Clust2MarkersMajor2+Clust2MarkersMinor2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.25, shape = 21,
              aes(fill = ifelse(InsufficientEvidence == 0, "Below cov thresholds", "Meets cov thresholds"))) +
  stat_compare_means() +
  geom_hline(yintercept = 7.5, linetype = 2) +
  xlab("Group 2 isolate obtained") +
  ylab("Number of Group 2 markers detected") + scale_fill_brewer(palette = "Set2", name = "")


pileup_results_summarize[,Clust2TM2:=ifelse(Clust2MarkersMajor+Clust2MarkersMinor==0, 0, Clust2TM)]
pileup_results_summarize[,Clust1TM2:=ifelse(Clust1MarkersMajor+Clust1MarkersMinor==0, 0, Clust1TM)]

pileup_results_summarize[InsufficientEvidence==0,table(Clust2TM2 > 4, NumClust2)]
pileup_results_summarize[InsufficientEvidence==0,table(Clust1TM2 > 4, NumClust1)]

# Compare with numbers in isolate genomes
genome_counts <- fread("snps_clust_markers_counts_per_genome.tsv")
genome_counts[Clust2==2,summary(Clust2Markers)]
genome_counts[Clust2==2, hist(Clust2Markers)]

# Save plots - Figure S6C
ggsave(clust2_markers_mgs, file = "mgs_method_im_dataset.pdf", width = 4, height = 3.5)

# Figure S6B
isolate_marker_plot2b <- ggplot(genome_counts[Clust2==2], aes(x = 1, y = HighPrevClust2Markers)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.25)+ xlab("Group B isolates") +
  ylab("Number of Group B markers") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(isolate_marker_plot2b, file = "groupB_snp_markers_isolates.pdf", width = 2.1, height = 2.7)
## Prevalence
priv_genes_all[ClustPrevalence > 0.1, table(ClusterMarker)]



library(pROC)
## Cutoffs for classifying Group B (clust2) present

preds <- pileup_results_summarize[InsufficientEvidence==0 & !is.na(NumClust2),Clust2MarkersMajor+Clust2MarkersMinor]
roc2 <- roc(response = pileup_results_summarize[InsufficientEvidence==0& !is.na(NumClust2), NumClust2 > 0],
            predictor = preds)
plot(roc2)
auc(roc2)
thresh2 <- sort(unique(pileup_results_summarize[InsufficientEvidence == 0 & !is.na(NumClust2),  Clust2MarkersMajor+Clust2MarkersMinor]))
thresh_tab <- data.table(thresh2, Accuracy = sapply(thresh2, function(x){
  pileup_results_summarize[InsufficientEvidence == 0 & !is.na(NumClust2)][,sum((Clust2MarkersMajor+Clust2MarkersMinor >= x) == (NumClust2 > 0))/length(NumClust2)]
}),
Sensitivity =   sapply(thresh2, function(x){
  pileup_results_summarize[InsufficientEvidence == 0 & !is.na(NumClust2) & NumClust2 > 0][,sum(Clust2MarkersMajor+Clust2MarkersMinor >= x)/length(NumClust2)]
}),
Specificity =   sapply(thresh2, function(x){
  pileup_results_summarize[InsufficientEvidence == 0 & !is.na(NumClust2) & NumClust2 == 0][,sum(Clust2MarkersMajor+Clust2MarkersMinor < x)/length(NumClust2)]
})
) ## Highest accuracy is at 17
# 8 allows for >93% sensntivity with 59% spec

pileup_results_summarize[,Clust1RelAbund:=Clust1Coverage*1e9/V2]
pileup_results_summarize[,Clust2RelAbund:=Clust2Coverage*1e9/V2]
pileup_results_summarize[,TotCovRelAbund:=TotalCovMean*1e9/V2]

## Compare with full Metaphlan data
pileup_results_summarize[,MetaphlanElenta:=NULL]
pileup_results_summarize <- merge(pileup_results_summarize, metaphlan_dat, by = "Sample", all.x = T)

# Figure S6D
metaphlan_compare <- ggplot(pileup_results_summarize, aes(x = MetaphlanElenta, y = TotCovRelAbund,
                                                          fill = ifelse(InsufficientEvidence == 0, "Below cov thresholds", "Meets cov thresholds"))) +
  geom_point(shape = 21, size = 2) + scale_fill_brewer(palette = "Set2", name = "") +
  xlab("E. lenta relative abundance estimated by Metaphlan4") + ylab("Read-normalized total marker coverage") +
  stat_cor()
ggsave(metaphlan_compare, file = "metaphlan_compare_markers.pdf", width = 5.6, height = 3.6)
pileup_results_summarize[InsufficientEvidence==0, cor.test(MetaphlanElenta, TotCovRelAbund,
                                                           method = "spearman")]
pileup_results_summarize[, cor.test(MetaphlanElenta, TotCovRelAbund,
                                    method = "spearman")]

## Final counts
pileup_results_summarize[InsufficientEvidence==0, table(Clust2MarkersMajor+Clust2MarkersMinor > 7, Clust1MarkersMajor+Clust1MarkersMinor > 7)]


#OK this is very consistent
pileup_results_summarize[,MetaphlanElentaCov:=as.numeric(MetaphlanElentaCov)]
coverage_compare_plot <- ggplot(pileup_results_summarize, aes(x = MetaphlanElentaCov, y = TotalCovMean)) + geom_point(shape = 21, alpha = 0.7, size = 1.7) + theme_classic2() +
  xlab("E. lenta read coverage (Metaphlan)") + ylab("Average coverage of all cluster markers") + stat_cor()


## Figure S6F
group2_total_compare_plot <- ggplot(pileup_results_summarize[InsufficientEvidence==0 & Clust2Present2==1], aes(x = Clust2Coverage, y = TotalCovMean)) + 
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray") +geom_point() + xlab("Group B coverage") +
  ylab("Total E. lenta coverage") + scale_x_sqrt() + scale_y_sqrt(breaks = c(10, 30, 50, 70, 100)) 
ggsave(group2_total_compare_plot, file = "immunomicrobiome_groupBcov_vs_total.pdf", width = 3.2, height = 3.2)

## Figure S6E - humann gene quantifications
### Graph Group 2 gene abundances vs Group 2 SNP abundances
gene_data <- fread("immunomicrobiomeMGS/genefamily_abundances_long.csv")
normalization_data <- fread("immunomicrobiomeMGS/normalization_data.txt")
gene_data <- merge(gene_data, normalization_data[,list(SampleID, ReadPairsPostQC, GenomeEquivalents)], by = "SampleID")

# Get key gene info from isolate analysis
top_gene_info <- fread("signifGenes_clust2_enrichmentImportanceSynteny2.tsv")
genes_of_interest <- fread("top_genes_annots_eggNOG_Uniref.tsv")
genes_of_interest_sub <- genes_of_interest[grepl("UniRef90", DbXrefs)]
genes_interest_1 <- fread("group1_assoc_genes_eggNOG.tsv")
genes_interest_core <- fread("core_genes_eggnog_Uniref.tsv")
genes_of_interest_sub[,Group:="Group 2 assoc"]
genes_interest_1[,Group:="Group 1 assoc"]
genes_interest_core[,Group:="Soft core"]
genes_of_interest_all <- rbind(genes_of_interest_sub, genes_interest_1, genes_interest_core, fill = T)
genes_of_interest_all <- distinct(genes_of_interest_all[,list(eggNOG_OGs, Uniref90, Group)])
gene_data[,GeneUniref90:=gsub("UniRef90_", "", Gene, fixed = T)]
gene_data <- gene_data[GeneUniref90 %in% genes_of_interest_all[,Uniref90]]
gene_list <- gene_data[,sort(unique(GeneUniref90))]
num_samps <- gene_data[,length(unique(SampleID)), by = GeneUniref90]
gene_data <- gene_data[GeneUniref90 %in% num_samps[V1 >=10, GeneUniref90]]
gene_list <- gene_data[,sort(unique(GeneUniref90))]

top_g2_genes <- top_gene_info[Importance > 0.004 | estimate > 4.5, GeneCluster]
gene_data_clust2 <- gene_data[GeneUniref90 %in% genes_of_interest_all[Group == "Group 2 assoc", Uniref90]]
gene_data_clust2 <- gene_data[GeneUniref90 %in% genes_of_interest_all[Group == "Group 2 assoc" & eggNOG_OGs %in% top_g2_genes, Uniref90]]

gene_list2 <- gene_data_clust2[,sort(unique(GeneUniref90))]
num_samps2 <- gene_data_clust2[,length(unique(SampleID)), by = GeneUniref90]

gene_data_clust2 <- dcast(gene_data_clust2, SampleID~GeneUniref90, value.var = "RPKG", fill = 0, fun.aggregate = sum) %>% 
  melt(id.var = "SampleID", variable.name = "GeneUniref90", variable.factor = F)

gene_data_clust2_compare <- merge(gene_data_clust2, pileup_results_summarize[,list(Sample, Clust2Coverage, Clust2RelAbund, 
                                                                                   InsufficientEvidence, Clust2Present, 
                                                                                   Clust2Present2, MetaphlanElenta, NumClust2)], 
                                  by.x = "SampleID", by.y = "Sample")

cor_results <- gene_data_clust2_compare[,
                                        list(cor(value, Clust2RelAbund, use = "complete.obs", method = "spearman"), cor.test(value, Clust2RelAbund, use = "complete.obs", method = "spearman")$p.value), by = list(GeneUniref90, GeneUniref90 %in% top_genes2)][order(V1, decreasing = T)]

cor_results <- merge(cor_results, distinct(genes_of_interest[,list(eggNOG_OGs, Description, Uniref90)]), by.x = "GeneUniref90", by.y = "Uniref90", all.x = T)
gene_data_clust2_compare <- merge(gene_data_clust2_compare, distinct(genes_of_interest[,list(eggNOG_OGs, Description, Uniref90)]), by.x = "GeneUniref90", by.y = "Uniref90", all.x = T)

tot_genes <- gene_data_clust2_compare[,sum(value), by = list(SampleID, InsufficientEvidence, Clust2Coverage, Clust2RelAbund, 
                                                             MetaphlanElenta, NumClust2)]
gene_snp_comparison <- ggplot(tot_genes, aes(y = Clust2RelAbund, x = V1, fill = factor(InsufficientEvidence))) + geom_point(shape = 21) + 
  scale_x_log10() + scale_fill_brewer(palette = "Set2", name = "") + stat_cor(method = "spearman") +
  xlab("Total RPKG of Group B-linked genes") + ylab("SNP-based relative abundance estimation")
ggsave(gene_snp_comparison, file = "immunomicrobiome_genes_snps_compare.pdf", width = 4.3, height = 3.6)

########### Figure S7 - compare with Cytof data

cytof <- fread("immuno_data/cytof.csv")
setnames(cytof, gsub("XHLT2-", "", names(cytof)))
cytof <- melt(cytof, id.var = "feature", variable.name = "Sample", variable.factor = F)
cytof <- merge(cytof, pileup_results_summarize, by = "Sample")
cytof[,FeatureCategory:=gsub("_[0-9]$", "", gsub("^[0-9]+\\.", "", feature))]
top_cats <- cytof[,length(unique(feature)), by = FeatureCategory][V1 > 6, FeatureCategory]
cytof[,FeatureCategory:=ifelse(FeatureCategory %in% top_cats, FeatureCategory, "Other")]
cytof[,list(
  cor(value, Clust2RelAbund, use = "complete.obs", method = "spearman"),
  cor.test(value, Clust2RelAbund, use = "complete.obs", method = "spearman")$p.value,
  cor(value, Clust1RelAbund, use = "complete.obs", method = "spearman"), 
  cor.test(value, Clust1RelAbund, use = "complete.obs", method = "spearman")$p.value, 
  cor(value, TotCovRelAbund, use = "complete.obs", method = "spearman"),
  cor.test(value, TotCovRelAbund, use = "complete.obs", method = "spearman")$p.value
), 
by = feature] %>% filter(V2 < 0.1|V4 < 0.1|V6 < 0.1)

cytof_assoc_results <- cytof[InsufficientEvidence==0,list(
  cor(value, Clust2RelAbund, use = "complete.obs", method = "spearman"),
  cor.test(value, Clust2RelAbund, use = "complete.obs", method = "spearman")$p.value,
  cor(value, Clust1RelAbund, use = "complete.obs", method = "spearman"), 
  cor.test(value, Clust1RelAbund, use = "complete.obs", method = "spearman")$p.value, 
  cor(value, TotCovRelAbund, use = "complete.obs", method = "spearman"),
  cor.test(value, TotCovRelAbund, use = "complete.obs", method = "spearman")$p.value
), 
by = list(feature, FeatureCategory)] #%>% filter(V2 < 0.1|V4 < 0.1|V6 < 0.1)
cytof_assoc_results[,Clust2PValAdj:=p.adjust(V2, method = "BH")]
cytof_assoc_results[,Clust1PValAdj:=p.adjust(V4, method = "BH")]
cytof_assoc_results[,TotElentaPValAdj:=p.adjust(V6, method = "BH")]

cytof_assoc_group2 <- ggplot(cytof_assoc_results, aes(x = V1, y = -1*log10(V2))) + 
  geom_hline(yintercept = -1*log10(0.05), linetype = 2) +
  geom_point(aes(fill = FeatureCategory), shape = 21, size = 2.5, alpha = 0.8) + xlim(-0.5, 0.5)+ 
  xlab("Spearman correlation") + ylab("-log10(p-value)") +
  ggtitle("CyTOF vs Group 2 E lenta\n(SNP-based)")+ 
  scale_fill_manual(values = c(brewer.pal(12, "Set3"), brewer.pal(3, "Dark2")), guide = guide_legend(ncol = 1))

cytof_assoc_tot <- ggplot(cytof_assoc_results, aes(x = V5, y = -1*log10(V6))) + 
  geom_hline(yintercept = -1*log10(0.05), linetype = 2) +
  geom_point(aes(fill = FeatureCategory), shape = 21, size = 2.5, alpha = 0.8) + xlim(-0.5, 0.5)+ 
  xlab("Spearman correlation") + ylab("-log10(p-value)") +
  ggtitle("CyTOF vs Total E lenta\n(SNP-based)")+ 
  scale_fill_manual(values = c(brewer.pal(12, "Set3"), brewer.pal(3, "Dark2")), guide = guide_legend(ncol = 1))


combined_plots2 <- cytof_assoc_tot+cytof_assoc_group2+ plot_layout(guides = "collect")
ggsave(combined_plots2, file = "cytof_plots_simple.pdf", width = 8.5, height = 3.5)

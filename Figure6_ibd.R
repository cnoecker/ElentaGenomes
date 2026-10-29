## Figure 6 - Metagenomic analysis - IBD data
library(data.table) 
library(tidyverse)
library(Biostrings)
library(readxl)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(ggrepel)
theme_set(theme_classic2())

metadata_qc <- readRDS("final_qc_metadata.rds")
el_metadata <- metadata_qc[grepl("lenta", species)]
priv_genes_all <- fread("snps_clust_markers_core_genes.tsv")

source("process_pileup_func.R")

### PRISM samples

# ################### PRISM data
samp_list <- list.files("prism_data/aligned_files/", pattern = "pileup_out.txt", full.names = T)
samp_list <- samp_list[grepl("SRR", samp_list)]
prism_pileup1 <- rbindlist(lapply(samp_list, process_pileup, cov_threshold = 0))
# load("mgs_processed_prism_covThresh0.rda")
# prism_pileup1 <- pileup_results

## HMBR data
# ################### Rest of IBD data
samp_list <- list.files("hmbr_data/aligned_files/", pattern = "pileup_out.txt", full.names = T)
pileup_results <- rbindlist(lapply(samp_list, process_pileup, cov_threshold = 0))
# load("mgs_processed_hmbrIBD_covThresh0.rda")

## Combine the 2 datasets
pileup_combined <- rbind(prism_pileup1, pileup_results, fill = T)
pileup_combined[,DownloadSource:=ifelse(grepl("SRR", Sample, fixed = T), "PRISM study/NCBI", "HMBR")]

## Combine sample metadata
# Old stuff first
supp_dat <- read_xlsx("prism_data/41564_2018_306_MOESM6_ESM.xlsx") %>% data.table()
names(supp_dat) <- c("Feature", gsub("|", "_", unlist(supp_dat[1, 2:ncol(supp_dat)]), fixed = T))
supp_dat_metadata <- supp_dat[1:9] %>% t() %>% data.table()
names(supp_dat_metadata) <- c("Feature", unlist(supp_dat_metadata[1,2:ncol(supp_dat_metadata)]))
supp_dat_metadata <- supp_dat_metadata[2:nrow(supp_dat_metadata)]
supp_dat_metadata[,Feature:=gsub("|", "_", Feature, fixed = T)]
sra_run_table <- fread("prism_data/SraRunTable.txt")
sra_run_table2 <- fread("prism_data/SraRunTableBigPRISM.txt")
sra_run_table <- rbind(sra_run_table[,list(Run, `Library Name`, `Sample Name`, AvgSpotLen, Bases)], 
                       sra_run_table2[,list(Run, `Library Name`, `Sample Name`, AvgSpotLen, Bases, AGE, sex, subject, diagnosis)], fill = T)
sra_run_table[,Source:=ifelse(is.na(AGE), "2018 study", "2024 study")]
sra_run_table <- merge(sra_run_table, supp_dat_metadata, by.x = "Library Name", by.y = "SRA_metagenome_name", all.x = T)
sra_run_table[,DownloadSource:="PRISM study/NCBI"]

# Kumbhari table S1
kumb_s1 <- read_xlsx("mmc2(1).xlsx") %>% data.table()
kumb_s1[,sample_id2:=gsub("[a-z0-9]+\\.", "", `Sample ID`)]
## New stuff - from HMBR
hmbr_metadata <- fread("hmbr_data/hmbr_compiled.tsv")
hmbr_metadata[,File_ID:=gsub("_1.fastq.gz", "", basename(R1))]
hmbr_metadata[,File_ID:=gsub(".adapterTrimmed.1_kneaddata_paired", "", File_ID)]
hmbr_metadata[,File_ID2:=ifelse(grepl("SKS", File_ID), gsub(".*_", "", File_ID), File_ID)]
hmbr_metadata[,table(File_ID2 %in% sra_run_table[,Run])]
hmbr_metadata[File_ID2 %in% sra_run_table[,Run], table(study_name)] ## CS-PRISM

hmbr_metadata[study_name=="LSS-PRISM", sample_accession:=gsub("_", "", paste0("p", sample_accession))]
# Data checks
hmbr_metadata[,table(sample_accession %in% kumb_s1[,sample_id2], study_name)] ## OK this looks good now
hmbr_metadata <- merge(hmbr_metadata, kumb_s1, by.x = "sample_accession", by.y = "sample_id2", all.x = T)
hmbr_metadata <- hmbr_metadata[,list(sample_accession, sample_id, study_name, subject_accession, gender, phenotype, race, bmi, 
                                     age_category, age_at_diagnosis_category, family_id, calprotectin, immunosuppressants, 
                                     steroids, antibiotics, database, `5ASA`, timepoint, File_ID, `Sample ID`, 
                                     `Study cohort`, Subject, Diagnosis, `Age (years)`, Sex, `BMI (kg/m^2)`,
                                     `Fecal calprotectin`, `Time relative to first sample (day; if longitudinal)`)]

sra_run_table_sub <- sra_run_table[!`Library Name` %in% hmbr_metadata[,sample_accession] & !Run %in% hmbr_metadata[,File_ID]]
overlapping_samps <- sra_run_table[`Library Name` %in% hmbr_metadata[,sample_accession] | Run %in% hmbr_metadata[,File_ID]]

sra_run_table_sub[is.na(diagnosis), diagnosis:=Diagnosis]
setnames(sra_run_table_sub, c("AGE", "sex", "subject", "diagnosis", "Fecal.Calprotectin", "antibiotic", "immunosuppressant"), 
         c("Age (years)", "Sex", "Subject", "phenotype", "calprotectin", "antibiotics", "immunosuppressants"))
metadata_all <- rbind(hmbr_metadata, sra_run_table_sub, fill = T)
metadata_all[is.na(DownloadSource), DownloadSource:="HMBR"]
metadata_all[is.na(sample_id), sample_id:=Run]

# Combine samples that appear twice
overlapping_samps[,IDMerge:=ifelse(`Library Name` %in% hmbr_metadata[,sample_accession], `Library Name`, Run)]
metadata_all[sample_accession %in% overlapping_samps[,`Library Name`]|File_ID %in% overlapping_samps[,Run],
             IDMerge:=ifelse(sample_accession %in% overlapping_samps[,`Library Name`], sample_accession, File_ID)]
metadata_all <- merge(metadata_all, overlapping_samps[IDMerge %in% metadata_all[,IDMerge], list(IDMerge, Bases, sex, subject, diagnosis, 
                                                                                                Age, Diagnosis, `Fecal.Calprotectin`, 
                                                                                                antibiotic, immunosuppressant, steroids)], 
                      by = "IDMerge", all.x = T)

metadata_all[,DiagnosisFinal:=case_when(
  (Diagnosis.x == "CD"|is.na(Diagnosis.x)) & phenotype=="CD" ~ "CD",
  (Diagnosis.x == "UC"|is.na(Diagnosis.x)) & phenotype=="UC" ~ "UC",
  (Diagnosis.x %in% c("Healthy", "Control")|is.na(Diagnosis.x)) & phenotype %in% c("HC", "Control", "Healthy") ~ "HC", 
  phenotype == "Indeterminate colitis" ~ "IC"
)]
metadata_all[,SexFinal:=case_when(
  Sex %in% c("F", "female") ~ "F",
  Sex %in% c("M", "male") ~ "M", 
  TRUE ~ "Missing"
)]
metadata_all[,AgeCategoryFinal:=case_when(
  !is.na(age_category) ~ age_category, 
  `Age (years)`  < 13 ~ "Child (<13 years)", 
  `Age (years)`  < 18 ~ "Adolescent (<18yrs)", 
  `Age (years)` < 65 ~ "Adult (<65yrs)",
  !is.na(`Age (years)`) ~ "Senior (65+)", 
  Age.x  < 13 ~ "Child (<13 years)", 
  Age.x  < 18 ~ "Adolescent (<18yrs)", 
  Age.x < 65 ~ "Adult (<65yrs)",
  !is.na(Age.x) ~ "Senior (65+)", 
  TRUE ~ "Missing"
)]
metadata_all[,antibiotics:=case_when(
  antibiotics %in% c("TRUE", "Yes") ~ "Yes",
  antibiotics %in% c("FALSE", "No") ~ "No", 
  TRUE ~ "Missing"
)]
metadata_all[,immunosuppressants:=case_when(
  immunosuppressants %in% c("True", "Yes") ~ "Yes",
  immunosuppressants %in% c("False", "No") ~ "No", 
  TRUE ~ "Missing"
)]
metadata_all[,calprotectinFinal:=case_when(
  !is.na(as.numeric(calprotectin)) ~ as.numeric(calprotectin), 
  !is.na(as.numeric(`Fecal calprotectin`)) ~ as.numeric(`Fecal calprotectin`), 
  !is.na(as.numeric(Fecal.Calprotectin)) ~ as.numeric(Fecal.Calprotectin)
)]
metadata_all[,subjectFinal:=ifelse(is.na(subject_accession), Subject, subject_accession)]
metadata_all[is.na(calprotectinFinal) & (sample_accession %in% sra_run_table[Fecal.Calprotectin != "None",`Library Name`]), 
             table(study_name)] #NL-IBD also duplicates, apparently - fixed
setnames(metadata_all, "Time relative to first sample (day; if longitudinal)", "TimeFromFirstSamp")

## Filter to just the samples with metadata
pileup_combined <- pileup_combined[Sample %in% metadata_all[,sample_id]]
pileup_combined[,length(unique(Sample))]

# Save SRA IDs for running Humann
sra_run_table_sub2 <- sra_run_table_sub[Run %in% pileup_combined[,Sample]]
fwrite(sra_run_table_sub2[,list(Run, `Library Name`, `Sample Name`, Bases, Source)], file = "prism_data/SraRunTableFinalSubset.txt", sep = ",")

## Summarize pileup data
pileup_combined[,Clust2Major:=ifelse(ClusterMarker == "Clust2" & MajorAllele == SNPMarker, 1, 0)]
pileup_combined[,Clust2Minor:=ifelse(ClusterMarker == "Clust2" & MinorAllele == SNPMarker & !is.na(MinorAllele), 1, 0)]
pileup_combined[,Clust1Major:=ifelse(ClusterMarker == "Clust1" & MajorAllele == SNPMarker, 1, 0)]
pileup_combined[,Clust1Minor:=ifelse(ClusterMarker == "Clust1" & MinorAllele == SNPMarker & !is.na(MinorAllele), 1, 0)]

pileup_combined[,Clust2Coverage:=case_when(
  Clust2Major == 1 ~ MajorAlleleCount,
  Clust2Minor == 1 ~ MinorAlleleCount,
  TRUE ~ 0
)]
pileup_combined[,Clust1Coverage:=case_when(
  Clust1Major == 1 ~ MajorAlleleCount,
  Clust1Minor == 1 ~ MinorAlleleCount,
  TRUE ~ 0
)]
pileup_results_summarize <- pileup_combined[,list(
  Clust2MarkersMajor = length(MajorAlleleCount[ClusterMarker == "Clust2" & MajorAlleleCount > 4 & MajorAllele == SNPMarker]),
  Clust2MarkersMinor = length(MinorAlleleCount[ClusterMarker == "Clust2" & MinorAlleleCount > 4 & MinorAllele == SNPMarker]),
  Clust2MarkersMajor2 = length(MajorAlleleCount[ClusterMarker == "Clust2" & MajorAlleleCount > 3 & MajorAllele == SNPMarker]),
  Clust2MarkersMinor2 = length(MinorAlleleCount[ClusterMarker == "Clust2" & MinorAlleleCount > 3 & MinorAllele == SNPMarker]),
  Clust2MedianCount = median(Clust2Coverage[Clust2Major==1 |Clust2Minor==1], na.rm = T),
  Clust2Q25Count = quantile(Clust2Coverage[Clust2Major==1 |Clust2Minor==1], na.rm = T, probs = 0.25),
  Clust2Q75Count = quantile(Clust2Coverage[Clust2Major==1 |Clust2Minor==1], na.rm = T, probs = 0.75),
  Clust2TM = mean(Clust2Coverage[Clust2Major == 1| Clust2Minor == 1], na.rm = T, trim = 0.25),
  Clust2Coverage = sum(Clust2Coverage[Clust2Major==1 |Clust2Minor==1])/length(Clust2Coverage[(Clust2Major==1 |Clust2Minor==1)]),
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
  TotGenesWCoverage2 = length(unique(Gene[NumCharSimplified > 2]))
), by = Sample]

pileup_results_summarize[,Clust1IQR:=Clust1Q75Count-Clust1Q25Count]
pileup_results_summarize[,Clust2IQR:=Clust2Q75Count-Clust2Q25Count]
pileup_results_summarize[,Clust1Present:=ifelse(Clust1MarkersMajor + Clust1MarkersMinor >= 8, 1, 0)]
pileup_results_summarize[,Clust2Present:=ifelse(Clust2MarkersMajor + Clust2MarkersMinor >= 8, 1, 0)]

pileup_results_summarize[,Clust1Present2:=ifelse(Clust1MarkersMajor2 + Clust1MarkersMinor2 >= 8, 1, 0)]
pileup_results_summarize[,Clust2Present2:=ifelse(Clust2MarkersMajor2 + Clust2MarkersMinor2 >= 8, 1, 0)]
pileup_results_summarize[,InsufficientEvidence:=ifelse(TotGenesWCoverage < 21, 1, 0)]

pileup_results_summarize[InsufficientEvidence == 0,Clust12Ratio:=case_when(
  Clust1Present == 1 & Clust2Present == 1 ~ Clust1TM/Clust2TM,
  Clust1Present == 1 & Clust2Present == 0 ~ 10,
  Clust1Present == 0 & Clust2Present == 1 ~ 0,
)]

## Combine with metadata
pileup_results_summarize <- merge(pileup_results_summarize, metadata_all[,list(sample_id, study_name, DiagnosisFinal, 
                                                                               SexFinal, AgeCategoryFinal, calprotectinFinal, 
                                                                               antibiotics, immunosuppressants, subjectFinal, 
                                                                               TimeFromFirstSamp, timepoint)], 
                                  by.x = "Sample", by.y = "sample_id")

pileup_results_summarize[,Outcome2:=case_when(
  TotGenesWCoverage2 < 21 ~ "InsufficientElentaEv", 
  Clust1Present2==1 & Clust2Present2 == 0 ~ "Clust1 only",
  Clust1Present2 == 0 & Clust2Present2 == 1 ~ "Clust2 only",
  Clust1Present2 == 1 & Clust2Present2 == 1 ~"Both clusters",
  Clust1Present2 == 0 & Clust2Present2 == 0 ~ "Neither/ambiguous"
)]

pileup_results_summarize[,Outcome2b:=case_when(
  Outcome2 == "InsufficientElentaEv" ~ "Insufficient E. lenta coverage", 
  Outcome2 %in% c("Both clusters", "Neither/ambiguous") ~ "Evidence of both/ambiguous", 
  Outcome2 == "Clust1 only" ~ "Group 1 present",
  Outcome2 == "Clust2 only" ~ "Group 2 present"
)]
pileup_results_summarize[,Outcome4:=case_when(
  Outcome2 == "InsufficientElentaEv" ~ "Insufficient E. lenta coverage",
  Clust2Present2==1 ~ "Lineage B present",
  TRUE ~ "Other E. lenta present (B not detected)"
)]

#Combine with read counts
read_counts_prism <- list.files("prism_data/humann/", pattern = "count.txt")
read_counts_hmbr <- list.files("hmbr_data/humann/", pattern = "count.txt")
read_counts_prism <- rbindlist(lapply(read_counts_prism, function(x){
  foo <- fread(paste0("prism_data/humann/", x))
  foo <- tail(foo, 1)
  foo[,Sample:=gsub("_readcount.txt", "", x)]
  return(foo)
}))
read_counts_hmbr <- rbindlist(lapply(read_counts_hmbr, function(x){
  foo <- fread(paste0("hmbr_data/humann/", x))
  foo <- tail(foo, 1)
  foo[,Sample:=gsub("_readcount.txt", "", x)]
  return(foo)
}))
read_counts <- rbind(read_counts_prism, read_counts_hmbr)
setnames(read_counts, "V1", "TotalReads")
pileup_results_summarize <- merge(pileup_results_summarize, read_counts, by = "Sample", all.x = T)

##### Figure 6A-B
color_scale3 <- c("gray80", brewer.pal(7, "Blues")[c(7, 4)]) #color_scale2[c(1,4,2)]
names(color_scale3) <- c("Insufficient E. lenta coverage", "Lineage B present", 
                         "Other E. lenta present (B not detected)")

## Normalized abundance estimates - include samples with 0 Group B
pileup_results_summarize[Outcome4 != "Insufficient E. lenta coverage", 
                         ElentaAbund:=TotalCovMean*1e9/TotalReads]
pileup_results_summarize[Outcome4 != "Insufficient E. lenta coverage", 
                         GroupBAbund:=case_when(
                           Clust2MarkersMajor + Clust2MarkersMinor == 0 ~ 0, 
                           TRUE ~ Clust2TM*1e9/TotalReads)]

## Add a total reads min threshold for everything - clearly skews abundance estimates

group_b_abunds <- ggplot(pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage" & TotalReads > 1.5e7], 
                         aes(x = factor(DiagnosisFinal, levels = rev(c("UC", "CD", "HC"))), y = GroupBAbund))+ geom_boxplot(outlier.shape = NA)  + 
  geom_jitter(height =0, width = 0.25, alpha = 0.3, shape = 21, fill = color_scale3[2]) + geom_pwc(label = "p.adj") +
  ylab("Group B normalized abundance") + scale_x_discrete(labels = c("Healthy\ncontrols", "Crohn's\ndisease", "Ulcerative\ncolitis"), name = "") +
  scale_y_sqrt()
tot_elenta_abunds <- ggplot(pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage" & TotalReads > 1.5e7], 
                            aes(x = factor(DiagnosisFinal, levels = rev(c("UC", "CD", "HC"))), y = ElentaAbund))+ geom_boxplot(outlier.shape = NA)  + 
  geom_jitter(height =0, width = 0.25, alpha = 0.3, shape = 21, fill = color_scale3[3]) + geom_pwc(label = "p.adj") +
  ylab("E. lenta normalized abundance")  + scale_x_discrete(labels = c("Healthy\ncontrols", "Crohn's\ndisease", "Ulcerative\ncolitis"), name = "") +
  scale_y_sqrt()

ggsave(tot_elenta_abunds+group_b_abunds, file = "ibd_snpBased_elenta_groupB_abunds.pdf", width = 7, height = 3.7)

group_b_abunds <- ggplot(pileup_results_summarize[DiagnosisFinal != "IC" & 
                                                    Outcome4 != "Insufficient E. lenta coverage"], 
                         aes(x = factor(DiagnosisFinal, levels = rev(c("UC", "CD", "HC"))), y = GroupBAbund))+ 
  geom_boxplot(outlier.shape = NA)  + 
  geom_jitter(height =0, width = 0.25, alpha = 0.3, shape = 21, fill = color_scale3[2]) + 
  geom_pwc(label = "p.adj") +
  ylab("Group B normalized abundance") + scale_x_discrete(labels = c("Healthy\ncontrols", "Crohn's\ndisease", "Ulcerative\ncolitis"), name = "") 
tot_elenta_abunds <- ggplot(pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage"], 
                            aes(x = factor(DiagnosisFinal, levels = rev(c("UC", "CD", "HC"))), y = ElentaAbund))+ geom_boxplot(outlier.shape = NA)  + 
  geom_jitter(height =0, width = 0.25, alpha = 0.3, shape = 21, fill = color_scale3[3]) + geom_pwc(label = "p.adj") +
  ylab("E. lenta normalized abundance")  + scale_x_discrete(labels = c("Healthy\ncontrols", "Crohn's\ndisease", "Ulcerative\ncolitis"), name = "") 

#ggsave(tot_elenta_abunds+group_b_abunds, file = "ibd_snpBased_elenta_groupB_abunds.pdf", width = 7, height = 3.7)

pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage" & TotalReads > 1.5e7, glm(DiagnosisFinal=="HC"~log1p(GroupBAbund), family = "binomial")] %>% broom::tidy()
pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage" & TotalReads > 1.5e7, glm(DiagnosisFinal=="HC"~log1p(ElentaAbund), family = "binomial")] %>% broom::tidy()

pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage", lm(log1p(GroupBAbund)~factor(DiagnosisFinal, levels = c("HC", "CD", "UC"))+TotalReads)] %>% broom::tidy()
pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage", lm(log1p(ElentaAbund)~factor(DiagnosisFinal, levels = c("HC", "CD", "UC"))+TotalReads)] %>% broom::tidy()

pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage" & TotalReads > 1.5e7, lm(log1p(GroupBAbund)~factor(DiagnosisFinal, levels = c("HC", "CD", "UC"))+SexFinal+AgeCategoryFinal+immunosuppressants)] %>% broom::tidy()
pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage" & TotalReads > 1.5e7, lm(log1p(ElentaAbund)~factor(DiagnosisFinal, levels = c("HC", "CD", "UC"))+SexFinal+AgeCategoryFinal+immunosuppressants)] %>% broom::tidy()

pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage" & TotalReads > 1.5e7, lm(log1p(GroupBAbund)~factor(DiagnosisFinal, levels = c("HC", "CD", "UC")))] %>% broom::tidy()
pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage" & TotalReads > 1.5e7, lm(log1p(ElentaAbund)~factor(DiagnosisFinal, levels = c("HC", "CD", "UC")))] %>% broom::tidy()

pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage" & TotalReads > 1.5e7, fisher.test((DiagnosisFinal=="HC"), Outcome4)]
## Mean/se lollipop plot with estimate effect sizes? ?

library(lmerTest)
mod1 <- pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage" & TotalReads > 1.5e7, 
                                     lmer(log1p(GroupBAbund)~factor(DiagnosisFinal, levels = c("HC", "CD", "UC"))+(1|subjectFinal))]
mod2 <- pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage" & TotalReads > 1.5e7, 
                                     lmer(log1p(ElentaAbund)~factor(DiagnosisFinal, levels = c("HC", "CD", "UC"))+(1|subjectFinal))] 

mod_tab <- rbind(data.frame(summary(mod1)$coefficients, Model = "Group B abundance"),
      data.frame(summary(mod2)$coefficients, Model = "Total E. lenta abundance")) %>% data.table()
mod_tab[,Term:=rep(c("Intercept", "CD", "UC"), 2)]

names(mod_tab)[1:2] <- c("estimate", "stderr")
names(mod_tab)[5] <- "pval"
lmer_coef_plot <- ggplot(mod_tab[Term != "Intercept"], 
                         aes(x = estimate, y = Term, fill = Model, group = Model)) + 
  geom_segment(aes(x = estimate - stderr, xend = estimate + stderr, color = Model), position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(width = 0.2), size = 3, shape = 21) + 
  geom_text(aes(x = 1.4, label = paste0("p = ", round(pval, 3))), color = "black", position = position_dodge(width = 0.2))+
  xlim(0, 1.5) + xlab("Effect size of association\nwith log taxon abundance") + 
  scale_fill_manual(values = brewer.pal(7, "Blues")[c(7, 4)], guide = guide_legend(ncol = 1)) + 
  scale_color_manual(values = brewer.pal(7, "Blues")[c(7, 4)], guide = guide_legend(ncol = 1)) #+
  theme(legend.position = "bottom")
ggsave(lmer_coef_plot, file = "group_b_totElen_abund_assoc_coef.pdf", width = 6, height = 3.8)

mod_tab[Term != "Intercept"]

pileup_results_summarize[Outcome4 != "Insufficient E. lenta coverage" & TotalReads > 1.5e7,table(DiagnosisFinal)]

pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage", 
                         glmer(factor(DiagnosisFinal=="HC") ~ log1p(GroupBAbund)+log1p(ElentaAbund)+(1|subjectFinal), family = "binomial")] %>% summary()

pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage", list(
  mean(ElentaAbund), 
  median(ElentaAbund) 
  ), by = DiagnosisFinal]
pileup_results_summarize[DiagnosisFinal != "IC" & Outcome4 != "Insufficient E. lenta coverage", list(
  mean(GroupBAbund), 
  median(GroupBAbund) 
), by = DiagnosisFinal]

### Figure 6E Calprotectin
calprotectin_plot_outcome4 <- ggplot(pileup_results_summarize[DiagnosisFinal != "HC" & TotalReads > 1.5e7], 
                                     aes(x = factor(Outcome4), y = as.numeric(calprotectinFinal))) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height = 0, width = 0.15, shape = 21, aes(fill = DiagnosisFinal), size = 2) +
  theme_classic2() + geom_pwc(method = "dunn_test", y.position = 550, step.increase = 0.04, bracket.nudge.y = 0, label = "p.adj")+ xlab("") + ylab("Fecal calprotectin (mg/g feces)") + 
  scale_fill_brewer(palette = "Set3", 
                    labels = c("Crohn's disease", "Indeterminate colitis", "Ulcerative colitis"), name = "", guide = guide_legend(ncol = 1)) + ylim(0, 800) + theme(legend.position = "bottom")
ggsave(calprotectin_plot_outcome4, file = "prism_outcome4_calprotectin_readThresh.pdf", width = 3.8, height = 5)

# Use these p-values which include the outlier - Dunn's test with multiple hypothesis correction
ggplot(pileup_results_summarize[DiagnosisFinal != "HC"& TotalReads > 1.5e7], 
       aes(x = factor(Outcome4), y = as.numeric(calprotectinFinal))) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height = 0, width = 0.15, shape = 21, aes(fill = DiagnosisFinal), size = 2) +
  theme_classic2() + geom_pwc(method = "dunn_test", y.position = 600, step.increase = 0.06, bracket.nudge.y = 0, label = "p.adj")+ xlab("") + ylab("Fecal calprotectin (mg/g feces)") + 
  scale_fill_brewer(palette = "Set3", 
                    labels = c("Crohn's disease", "Indeterminate colitis", "Ulcerative colitis"), name = "", guide = guide_legend(ncol = 1)) + theme(legend.position = "bottom")

## Effect sizes
pileup_results_summarize[DiagnosisFinal != "HC"& TotalReads > 1.5e7, 
                         median(calprotectinFinal, na.rm = T), by = factor(Outcome4)]

pileup_results_summarize[,study2:=case_when(
  grepl("PRISM", study_name) ~ "PRISM", 
  grepl("NL", study_name) ~ "NL-IBD", 
  is.na(study_name) ~ "Other",
  TRUE ~ "Other"
)]
pileup_results_summarize[DiagnosisFinal != "HC" & DiagnosisFinal != "IC" & InsufficientEvidence==0 & TotalReads > 1.5e7, 
                         broom::tidy(lm(calprotectinFinal~factor(Clust2Present2)+factor(Clust1Present2)+DiagnosisFinal+
                                          study2))]

#### Humann analysis of associated genes - Figure 6 CD and FG
genes_of_interest <- fread("top_genes_annots_eggNOG_Uniref.tsv")
genes_of_interest_sub <- genes_of_interest[grepl("UniRef90", DbXrefs)]
genes_interest_1 <- fread("group1_assoc_genes_eggNOG.tsv")
genes_interest_core <- fread("core_genes_eggnog_Uniref.tsv")
genes_of_interest_sub[,Group:="Group B assoc"]
genes_interest_1[,Group:="Group A assoc"]
genes_interest_core[,Group:="Soft core"]
genes_of_interest_all <- rbind(genes_of_interest_sub, genes_interest_1, genes_interest_core, fill = T)
genes_of_interest_all <- distinct(genes_of_interest_all[,list(eggNOG_OGs, Uniref90, Group)])

## Elenmrk1 = FKHGEB_06550 in DSM 2243
elnmrk_uniref90 <- "F0HLY6"
genes_of_interest_all[Uniref90==elnmrk_uniref90]

# Import humann data
humann_files <- c(list.files("prism_data/humann/", pattern = "genefamilies.tsv", full.names = T), 
                  list.files("hmbr_data/humann/", pattern = "genefamilies.tsv", full.names = T))
humann_files <- humann_files[gsub("_genefamilies.tsv", "", basename(humann_files)) %in% metadata_all[,sample_id]]
humann_files <- humann_files[gsub("_genefamilies.tsv", "", basename(humann_files)) %in% read_counts[TotalReads > 1.5e7, Sample]]
humann_dat <- rbindlist(lapply(humann_files, function(x){
  foo <- fread(x)
  if(nrow(foo) > 0){
    foo[,c("GeneFamily", "Taxon"):=tstrsplit(`# Gene Family`, split = "\\|")]
    foo[,GeneFamily:=gsub("UniRef90_", "", GeneFamily)]
    setnames(foo, names(foo)[2],  "AbundanceRPKs")
    foo_interest <- foo[GeneFamily %in% genes_of_interest_all[,Uniref90]]
    foo_interest[,Sample:=gsub("_genefamilies.tsv", "", basename(x))]
    return(foo_interest)
  }
}))
humann_dat <- merge(humann_dat, read_counts, by = "Sample", all.x = T)
humann_dat[,RPKM:=AbundanceRPKs*1e6/TotalReads]
saveRDS(humann_dat, file = "humann_mgs_dat.rds")

humann_dat <- readRDS("humann_mgs_dat.rds")
genes_interest_sub <- distinct(genes_of_interest_all[,list(Uniref90, eggNOG_OGs, Group)]) %>% 
  filter(!Uniref90 %in% genes_of_interest_all[,length(unique(eggNOG_OGs)), by = Uniref90][V1 > 1, Uniref90]) ## Remove dups

humann_dat <- merge(humann_dat, genes_interest_sub, by.x = "GeneFamily", by.y = "Uniref90", 
                    all.x = T)

humann_dat <- humann_dat[TotalReads != 0] # Remove problematic sample
# Stratified data only
humann_dat <- humann_dat[!is.na(Taxon)]
humann_dat_zeros <- humann_dat %>% dcast(GeneFamily+Taxon+eggNOG_OGs+Group~Sample, value.var = "RPKM", fill = 0) %>% 
  melt(id.var = c("GeneFamily", "Taxon", "eggNOG_OGs", "Group"), variable.name = "Sample", variable.factor = F)

humann_dat_zeros <- merge(humann_dat_zeros, metadata_all[,list(sample_id, study_name, DiagnosisFinal, 
                                                               SexFinal, AgeCategoryFinal, calprotectinFinal, 
                                                               antibiotics, immunosuppressants, subjectFinal, 
                                                               TimeFromFirstSamp, timepoint)], 
                          by.x = "Sample", by.y = "sample_id")

# Make sure samples with 0 are included
# Apply same read count filter
humann_dat_zeros <- humann_dat_zeros[Sample %in% read_counts[TotalReads > 1.5e7, Sample]]
bad_genes <- humann_dat_zeros[!is.na(Group), sum(value), by = list(Sample, 
                                                                   eggNOG_OGs, GeneFamily, DiagnosisFinal, Group)
][,length(V1[V1 != 0]),by = list(eggNOG_OGs, GeneFamily, Group)][V1 < 4] # Present in at least 4 samples
humann_dat_zeros <- humann_dat_zeros[!GeneFamily %in% bad_genes[,GeneFamily]]

humann_model_results <- humann_dat_zeros[!is.na(Group), sum(value), by = list(Sample, study_name, eggNOG_OGs,GeneFamily, DiagnosisFinal, Group)] %>% 
  filter(Group != "Group 1 assoc") %>% 
  group_by(eggNOG_OGs, GeneFamily, Group) %>% 
  mutate(ScaledValue = scale(log(V1 + 0.0001)), 
         Disease = ifelse(DiagnosisFinal != "HC", 1, 0)) %>% 
  reframe(
    broom::tidy(glm(factor(Disease) ~ ScaledValue+study_name, family = "binomial"))
  ) %>% filter(grepl("ScaledValue", term)) %>% 
  mutate(Padj = p.adjust(p.value, method = "BH")) %>% 
  data.table() 


bad_genes_uc <- humann_dat_zeros[!is.na(Group) & DiagnosisFinal != "CD", sum(value), by = list(Sample, 
                                                                                               eggNOG_OGs, GeneFamily, DiagnosisFinal, Group)
][,length(V1[V1 != 0]),by = list(eggNOG_OGs, GeneFamily, Group)][V1 < 4] # Present in at least 4 samples

uc_model <- humann_dat_zeros[!is.na(Group), sum(value), by = list(Sample, study_name, eggNOG_OGs, GeneFamily, DiagnosisFinal, Group)] %>% 
  filter(Group != "Group 1 assoc") %>% 
  filter(!GeneFamily %in% bad_genes_uc[,GeneFamily]) %>% 
  group_by(eggNOG_OGs, GeneFamily, Group) %>% 
  filter(DiagnosisFinal %in% c("UC", "HC")) %>% 
  mutate(ScaledValue = scale(log(V1 + 0.0001)), 
         Disease = ifelse(DiagnosisFinal != "HC", 1, 0)) %>% 
  reframe(
    broom::tidy(glm(factor(Disease) ~ ScaledValue + study_name, family = "binomial"))
  ) %>% filter(grepl("ScaledValue", term)) %>% 
  mutate(Padj = p.adjust(p.value, method = "BH")) %>% 
  data.table() #%>% ggplot( aes(x = estimate, color = Group)) + geom_density()

bad_genes_cd <- humann_dat_zeros[!is.na(Group) & DiagnosisFinal != "UC", sum(value), by = list(Sample, 
                                                                                               eggNOG_OGs, GeneFamily, DiagnosisFinal, Group)
][,length(V1[V1 != 0]),by = list(eggNOG_OGs, GeneFamily, Group)][V1 < 4] # Present in at least 4 samples

cd_model <- humann_dat_zeros[!is.na(Group), sum(value), by = list(Sample, study_name, eggNOG_OGs, GeneFamily, DiagnosisFinal, Group)] %>% 
  filter(Group != "Group 1 assoc") %>% 
  filter(!GeneFamily %in% bad_genes_cd[,GeneFamily]) %>% 
  group_by(eggNOG_OGs, GeneFamily, Group) %>% 
  filter(DiagnosisFinal %in% c("CD", "HC")) %>% 
  mutate(ScaledValue = scale(log(V1 + 0.0001)), 
         Disease = ifelse(DiagnosisFinal != "HC", 1, 0)) %>% 
  reframe(
    broom::tidy(glm(factor(Disease) ~ ScaledValue + study_name, family = "binomial"))
  ) %>% filter(grepl("ScaledValue", term)) %>% 
  mutate(Padj = p.adjust(p.value, method = "BH")) %>% 
  data.table() 

top_gene_info <- fread("signifGenes_clust2_enrichmentImportanceSynteny2.tsv")
top_g2_genes <- top_gene_info[Importance > 0.004 | estimate > 4.5, GeneCluster]

uc_model <- merge(uc_model, top_gene_info[,list(GeneCluster, Description)], by.x = "eggNOG_OGs", by.y = "GeneCluster", all.x = T)
uc_model[,PointLab:=case_when(
  Description %in% c("HindVP restriction endonuclease", "DEAD-like helicases superfamily", "Mrr N-terminal domain",
                     "Cytosine-specific methyltransferase")|eggNOG_OGs=="2F3QR@1|root,33WHG@2|Bacteria,2H8D3@201174|Actinobacteria" ~ "Defense",
  Description %in% c("binding domain protein", "Nitroreductase family", "Bacterial regulatory proteins, tetR family", 
                     "Putative ATP-dependent DNA helicase recG C-terminal", "SpoVT / AbrB like domain") ~ "NR",
  TRUE ~ ""
)]

uc_assoc <- ggplot(uc_model, aes(x = estimate, y = -1*log10(p.value), fill = Group)) + 
  geom_hline(yintercept = 2, linetype = 2, color = "gray") +
  geom_point(data = uc_model[Group == "Soft core"], alpha = 0.1, shape = 21, color = "gray60", size = 1.6) +
  geom_point(data = uc_model[eggNOG_OGs %in% top_g2_genes], shape = 21, color = "black", size = 1.6) +
  scale_fill_manual(values = c(brewer.pal(7, "Blues")[7], "gray80"), 
                    labels = c("Group 2-linked genes", "E. lenta core genes"), name = "", guide = guide_legend(ncol = 1)) + 
  ylim(0, 8) + xlim(-1.05, 1.3) + ggtitle("Ulcerative colitis vs. healthy controls") +
  xlab("Log scaled association with disease") #+ 

cd_model <- merge(cd_model, top_gene_info[,list(GeneCluster, Description)], by.x = "eggNOG_OGs", by.y = "GeneCluster", all.x = T)

cd_model[,PointLab:=case_when(
  Description %in% c("HindVP restriction endonuclease", "DEAD-like helicases superfamily", "Mrr N-terminal domain",
                     "Cytosine-specific methyltransferase")|eggNOG_OGs=="2F3QR@1|root,33WHG@2|Bacteria,2H8D3@201174|Actinobacteria" ~ "Defense",
  Description %in% c("binding domain protein", "Nitroreductase family", "Bacterial regulatory proteins, tetR family", 
                     "Putative ATP-dependent DNA helicase recG C-terminal", "SpoVT / AbrB like domain") ~ "NR",
  TRUE ~ ""
)]
cd_assoc <- ggplot(cd_model, aes(x = estimate, y = -1*log10(p.value), fill = Group)) + 
  geom_hline(yintercept = 2, linetype = 2, color = "gray") +
  geom_point(data = cd_model[Group == "Soft core"], alpha = 0.1, shape = 21, color = "gray60", size = 1.6) +
  geom_point(data = cd_model[eggNOG_OGs %in% top_g2_genes], shape = 21, color = "black", size = 1.6) +
  scale_fill_manual(values = c(brewer.pal(7, "Blues")[7], "gray80"), 
                    labels = c("Group 2-linked genes", "E. lenta core genes"), name = "", guide = guide_legend(ncol = 1)) + 
  ylim(0, 8) + xlim(-1.25, 1.25)+ ggtitle("Crohn's disease vs. healthy controls") +  xlab("Log scaled association with disease") #+
#geom_text_repel(data = cd_model[PointLab != ""], aes(label = PointLab),min.segment.length = 0, size = 2.5, max.overlaps = 20, color = brewer.pal(7, "Blues")[7])

assoc_plot <- uc_assoc+cd_assoc + plot_layout(guides = "collect") #&
ggsave(assoc_plot, file = "mgs_ibd_humann_genes_lineageB_assoc_updated.pdf", width = 8.7, height = 3.4)

humann_dat_zeros <- merge(humann_dat_zeros, unique(top_gene_info[,list(GeneCluster, Description)]), by.x = "eggNOG_OGs", by.y = "GeneCluster", all.x = T)

ggplot(cd_model[Group != "Group A assoc"], aes(x = estimate, color = Group)) + geom_density()
ggplot(cd_model[Group != "Group A assoc"], aes(x = estimate, color = ifelse(eggNOG_OGs %in% top_g2_genes,T, F))) + geom_density()
ggplot(uc_model[Group != "Group A assoc"], aes(x = estimate, color = ifelse(eggNOG_OGs %in% top_g2_genes,T, F))) + geom_density()
ggplot(uc_model[Group != "Group A assoc"], aes(x = estimate, color = Group)) + geom_density()

## Correlations with calprotectin
corrs_all <- humann_dat_zeros[!is.na(Group), sum(value), 
                              by = list(Sample, study_name, eggNOG_OGs, GeneFamily, 
                                        DiagnosisFinal, Group, calprotectinFinal)] %>% 
  filter(Group != "Group 1 assoc" & !is.na(Group)) %>% 
  filter(!is.na(eggNOG_OGs)) %>% 
  group_by(eggNOG_OGs, Group) %>% 
  reframe(
    estimate = cor.test(V1,calprotectinFinal, method = "spearman")$estimate
  ) %>% 
  mutate(GeneGroup = ifelse(eggNOG_OGs %in% top_g2_genes, "Top Group B genes", "E. lenta core genes"))
corr_plot <- ggplot(corrs_all, aes(x = estimate, color = GeneGroup)) + geom_density(linewidth = 0.8) + 
  xlab("Spearman correlation with fecal calprotectin") + 
  scale_color_manual(values = c(brewer.pal(7, "Blues")[c(4, 7)]), name = "") + ylab("Density")
ggsave(corr_plot, file = "calprotectin_corr_dist_plot.pdf", width = 4.8, height = 3.2)  
#%>% ggplot( aes(x = estimate)) + geom_histogram()

## Correlations for all top genes
all_top_genes <- ggplot(humann_dat_zeros[eggNOG_OGs %in% top_g2_genes & DiagnosisFinal != "HC", sum(value), 
                                         by = list(Sample, eggNOG_OGs, calprotectinFinal, DiagnosisFinal, Description)], 
                        aes(y = calprotectinFinal, x = V1)) + geom_point(aes(fill = DiagnosisFinal), shape= 21) + 
  facet_wrap(~paste0(Description, "\n", eggNOG_OGs), scales = "free")+
  stat_cor(method = "spearman") + ylim(0, 750) + scale_fill_brewer(palette = "Set3", labels = c("Crohn's disease", "Indeterminate colitis", "Ulcerative colitis")) +
  scale_x_sqrt()

# Figure 6F-G - 2 representative Group B genes
top_genes_rep <- ggplot(humann_dat_zeros[(grepl("Nitroreduc", Description)|grepl("HindVP", Description))& DiagnosisFinal != "HC", sum(value), 
                                         by = list(Sample, eggNOG_OGs, calprotectinFinal, DiagnosisFinal, Description)], 
                        aes(y = calprotectinFinal, x = V1)) + 
  stat_smooth(method = "lm", color = "gray30", linewidth = 0.7, fill = "gray80", linetype = 2) +
  geom_point(aes(fill = DiagnosisFinal), shape= 21) + 
  facet_wrap(~Description, scales = "free")+ xlab("Gene abundance (RPKM)") +
  ylab("Fecal calprotectin (mg/g feces)") +
  stat_cor(method = "spearman") + ylim(0, 680) + 
  scale_fill_brewer(palette = "Set3", labels = c("Crohn's disease", "Indeterminate colitis", "Ulcerative colitis"), guide = "none") +
  scale_x_sqrt() + 
  theme(strip.background = element_blank())
ggsave(all_top_genes, file= "prism_toplineageBGenes_calprotectin.pdf", width = 10, height = 9)
ggsave(top_genes_rep, file= "prism_toplineageBGenesRep_calprotectin.pdf", width = 5.5, height = 3.5)

# Use this version to get correct p-values/corrs - ggpubr does not include off-scale data
ggplot(humann_dat_zeros[(grepl("Nitroreduc", Description)|grepl("HindVP", Description))& DiagnosisFinal != "HC", sum(value), 
                        by = list(Sample, eggNOG_OGs, calprotectinFinal, DiagnosisFinal, Description)], 
       aes(y = calprotectinFinal, x = V1)) + 
  stat_smooth(method = "lm", color = "gray30", alpha = 0.4) +
  geom_point(aes(fill = DiagnosisFinal), shape= 21) + 
  facet_wrap(~paste0(Description, "\n", eggNOG_OGs), scales = "free")+ xlab("Gene abundance (RPKM)") +
  ylab("Fecal calprotectin (mg/g feces)") +
  stat_cor(method = "spearman") + 
  scale_fill_brewer(palette = "Set3", labels = c("Crohn's disease", "Indeterminate colitis", "Ulcerative colitis")) +
  scale_x_sqrt()


## TNF data - import and summarize
el_metadata_full <- read_xlsx("Noecker_SupplTables_v1.xlsx", sheet ="TableS3", skip =2)

setnames(el_metadata_full, "Genome ID", "genome")
tnf_data_raw <- read_xlsx("sorted_ELISA_Results_26_07_23_RR.xlsx")

tnf_data_raw <- tnf_data_raw %>% mutate(
  Replicate = 1:nrow(tnf_data_raw)
) %>% 
  pivot_longer(-Replicate) %>% 
  filter(!is.na(value))

tnf_data_raw <- tnf_data_raw %>% 
  mutate(genome = case_when(
    grepl(10906, name) ~ "GCF_000185625.1_ASM18562v1_genomic",
    grepl(110905, name) ~  "GCF_003340245.1_ASM334024v1_genomic",
    grepl("DSM11767", name) ~ "GCF_003340045.1_ASM334004v1_genomic",
    grepl(2243, name) ~ "GCF_000024265.1_ASM2426v1_genomic",
    grepl("HS111", name) ~ "HS111-1A_CN11",
    grepl("HS44-R", name) ~ "HS44-R24B",
    TRUE ~ name
  ))

tnf_data_raw <- left_join(tnf_data_raw, el_metadata_full, by ="genome")

tnf_data_raw <- tnf_data_raw %>% 
  mutate(name  = factor(name, levels =c("PBS", sort(unique(tnf_data_raw$name))[sort(unique(tnf_data_raw$name))!="PBS"])))

tnf_data_summ <- tnf_data_raw %>% group_by(name, genome, `Genome group`) %>% 
  summarize(meanVal = mean(value), sdVal = sd(value), seVal = sd(value)/sqrt(length(value)))

strain_order <- tnf_data_summ %>% arrange(desc(meanVal)) %>% pull(name)

tnf_data_raw <- mutate(tnf_data_raw, 
                       `Genome group` = factor(ifelse(is.na(`Genome group`), "PBS", `Genome group`), levels = c("PBS", "A", "B")))

# Fig 6 H-I
plot_all_strains <- ggplot(tnf_data_raw, aes(x = factor(name, levels = strain_order), y = value, 
                         color =`Genome group`)) +
  geom_boxplot() +geom_jitter(height =0, width = 0.2, size = 0.5) +
  ylab("Cell count-normalized TNF secretion") + 
  theme(axis.text.x =element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_color_manual(values =c("gray", brewer.pal(7, "Blues")[c(4,7)])) +
  facet_grid(.~`Genome group`, space = "free", scales ="free_x") + 
  theme(strip.background = element_blank(), strip.text = element_text(size = 12), 
        legend.position = "bottom") +
  xlab("")
ggsave(plot_all_strains, file = "tnf_plot_all_strains.pdf", width = 6.7, height = 3.7)

tnf_plot <- ggplot(tnf_data_summ, aes(x = `Genome group`, y=meanVal)) + 
  geom_violin(draw_quantiles = 0.5, color = "gray30")+
  geom_pointrange(aes(fill = `Genome group`, ymin = meanVal-seVal, ymax = meanVal+seVal), shape = 21, 
                  position = position_jitter(height = 0, width = 0.2))+
  scale_fill_manual(values =c(brewer.pal(3, "Blues")[c(2,3)], "gray"), labels = c("A", "B", "PBS control")) +
  stat_compare_means(method = "wilcox", data = tnf_data_summ %>% filter(name != "PBS"), label = "p.format") +
  ylab("Cell count-normalized TNF secretion") + scale_x_discrete(labels = c("Group A", "Group B", "PBS control"), 
                                                                 name = "")

tnf_plot2 <- tnf_data_summ %>% filter(name != "PBS") %>% 
  ggplot(aes(x = `Genome group`, y=meanVal)) + 
  geom_violin(draw_quantiles = 0.5, color = "gray30")+
  geom_point(aes(fill = `Genome group`, y = meanVal), shape = 21, 
                  position = position_jitter(height = 0, width = 0.2), size = 2)+
  scale_fill_manual(values =c(brewer.pal(7, "Blues")[c(4,7)], "gray"), labels = c("A", "B")) +
  stat_compare_means(method = "wilcox", data = tnf_data_summ %>% filter(name != "PBS"), label = "p.format") +
  ylab("Cell count-normalized TNF secretion") + 
  scale_x_discrete(labels = c("Group A", "Group B"),name = "")


ggsave(tnf_plot2, file = "tnf_data_replot2.pdf", width = 4.2, height = 3.8)



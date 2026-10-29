## Figure 5 and Figure S5 - plasmid analysis
library(data.table)
library(tidyverse)
library(ggpubr)
library(Biostrings)
library(DECIPHER)
library(ComplexHeatmap)
library(RColorBrewer)
library(ape)
library(plasmapR)
library(patchwork)
library(ggrepel)
library(broom)
library(readxl)
library(vegan)
library(ggtree)

theme_set(ggpubr::theme_classic2())

metadata <- readRDS("final_qc_metadata.rds")
el_metadata <- metadata[grepl("lenta", species)]

## Filtering for confident plasmids - Genomad, Circle_map

genomad_results <- list.files("annots_out/genomad/")
genomad_results <- genomad_results[genomad_results %in% el_metadata[,genome]]
g_results <- rbindlist(lapply(genomad_results, function(x){
  foo <- fread(paste0("annots_out/genomad/", x, "/", x, "_aggregated_classification/", x, 
                      "_aggregated_classification.tsv"))
  foo[,genome:=x]
  return(foo)
}))
contig_lengths <- fread("all_elenta_contigs_lengths_bothIdsCorrect.tsv")

g_results <- merge(g_results, contig_lengths, by.x = c("genome", "seq_name"), by.y = c("genome", "contig"))
potential_plasmids <- g_results[plasmid_score > 0.5]
potential_plasmids <- merge(potential_plasmids, el_metadata[,list(genome, Clust2, IsMAG, Country, Continent, birthCountry, birthContinent, GenomeSource.x, collection_date, age, age_category, gender)])

# Combine different plasmid/contig names
plasmid_files2 <- list.files("annots_out/mash_dists_genomad/plasmids/")
plasmids_len2 <- rbindlist(lapply(plasmid_files2[3:length(plasmid_files2)], function(x){
  gen <- readDNAStringSet(paste0("annots_out/mash_dists_genomad/plasmids/", x))
  contig_tab <- data.table(genome = gsub("_[0-9]+.fna", "", x), plasmidname = gsub(".fna$", "", x), contig = names(gen), len = nchar(gen))
  return(contig_tab)
}), fill = T)
plasmids_len2 <- plasmids_len2[genome %in% el_metadata[, genome]]
plasmids_len2[,contigAlt:=sapply(contig, function(y){
  strsplit(y, split = " ")[[1]][2]
})]
setnames(plasmids_len2, c("contig", "contigAlt"), c("contigFull", "contig"))
plasmids_len2 <- merge(plasmids_len2, contig_lengths, by = c("genome", "len", "contig"))

## Circle map output
circle_map_files <- list.files("circleMapOut/")
shortreadfiles <- fread("hybrid_genomes_all.txt", header = F)[,gsub(".*\\/", "", V1)]
el_metadata[GenomeSource.x == "ImmunoMicrobiome" & !genome %in% gsub("_circle.bed", "", circle_map_files) &
              genome %in% shortreadfiles, genome]

circle_output <- rbindlist(lapply(circle_map_files, function(x){
  foo <- fread(paste0("circleMapOut/", x))
  setnames(foo, c("contig", "start", "end", "discordantReads", "splitReads", "circleScore", "meanCov", "stdevCov", 
                  "covIncreaseStart", "covIncreaseEnd", "covContinuity"))
  foo[,genome:=gsub("_circle.bed", "", x)]
  return(foo)
}))
# Only keeping the ones with metadata for now
circle_output <- merge(circle_output, contig_lengths, by = c("genome", "contig"))
circle_output[,circleSize:=end+1-start]
potential_plasmids2 <- circle_output[circleSize == len & circleScore > 100 & circleSize < 2e5]

## Get gene annots from bakta
bakta_files <- paste0(el_metadata[,genome], ".tsv")
bakta_dat <- rbindlist(lapply(bakta_files, function(x){
  foo <- fread(paste0("annots_out/bakta/", x))
  foo[,genome:=gsub(".tsv", "", x)]
  return(foo)
}))
setnames(bakta_dat, "#Sequence Id", "Contig")

## Define set of confidnet plasmids
potential_plasmids2 <- merge(potential_plasmids2, contig_lengths[,list(genome, contig, baktaContig)], by = c("genome", "contig"), all.x = T)
potential_plasmids2 <- merge(potential_plasmids2, select(g_results, -len), by.x = c("genome", "contig"), by.y = c("genome", "seq_name"), all.x = T)
potential_plasmids2 <- potential_plasmids2[chromosome_score < 0.9]
potential_plasmids2[,PlasmidID:=paste0(genome, baktaContig)]

bakta_im <- bakta_dat[genome %in% el_metadata[GenomeSource.x %in% c("ImmunoMicrobiome", "Harris Wang lab"), genome]]
bakta_im_plasmids <- bakta_im[paste0(genome, Contig) %in% potential_plasmids2[len < 100000,paste0(genome, baktaContig)]]

## Get plasmid genes, write to file
plasmid_genomes <- bakta_im_plasmids[,sort(unique(genome))]
dir.create("confident_plasmid_faas/")
all_plasmid_faas <- AAStringSet()
for(j in 1:length(plasmid_genomes)){
  faas <- readAAStringSet(paste0("finalset_faas/", plasmid_genomes[j], ".faa"))
  names(faas) <- gsub(" .*", "", names(faas))
  faas_sub <- faas[names(faas) %in% bakta_im_plasmids[,`Locus Tag`]]
  all_plasmid_faas <- append(all_plasmid_faas, faas_sub)
}
writeXStringSet(all_plasmid_faas, file = "all_confidentPlasmid_cds.faa")

bakta_im_plasmids[,PlasmidID:=paste0(genome, Contig)]
plasmid_set <- bakta_im_plasmids[,sort(unique(PlasmidID))]

for(j in 1:length(plasmid_set)){
  faas_sub <- all_plasmid_faas[names(all_plasmid_faas) %in% bakta_im_plasmids[PlasmidID==plasmid_set[j], `Locus Tag`]]
  writeXStringSet(faas_sub, file = paste0("confident_plasmid_faas/", plasmid_set[j], ".faa"))
}

el_eggnog <- readRDS("elenta_eggnog.rds")
top_cat <- el_eggnog[,length(unique(eggNOG_OGs)), by=COG_category][order(V1, decreasing = T)][1:10, COG_category]
el_eggnog[,COG_category2:=ifelse(COG_category %in% top_cat, COG_category, "Multiple/Other")]

plasmid_eggnog <- el_eggnog[`#query` %in% bakta_im_plasmids[,`Locus Tag`]]
plasmid_eggnog <- merge(plasmid_eggnog, bakta_im_plasmids[,list(Contig, Type, Start, Stop, Strand, `Locus Tag`, Gene, Product, PlasmidID)], by.x = "#query", by.y = "Locus Tag")
plasmid_eggnog <- merge(plasmid_eggnog, el_metadata[,list(genome, Clust2, GenomeSource.x, IsMAG, birthContinentMerged, Continent)], by = "genome")
plasmid_eggnog <- merge(plasmid_eggnog, contig_lengths[,list(genome, baktaContig, len)], by.x = c("genome", "Contig"), by.y = c("genome", "baktaContig"), all.x = T)

##Proteinortho-based clustering - group plasmids with similar coding sequences
po_plasmid <- fread("proteinortho_results/plasmidGenesElenta.proteinortho.tsv")
po_plasmid <- mutate(po_plasmid, ClusterID:=1:nrow(po_plasmid))
setnames(po_plasmid, gsub(".faa$", "", names(po_plasmid)))
po_plasmid <- melt(po_plasmid[,4:ncol(po_plasmid), with = F], id.var = "ClusterID", 
                   variable.name = "PlasmidID", variable.factor = F)
po_plasmid <- mutate(po_plasmid, Present = ifelse(value == "*", 0, 1))
pe_wide3 <- po_plasmid %>% dcast(PlasmidID~ClusterID, value.var = "Present", fill = 0, fun.aggregate = sum)
plasmids_gene_dist <- vegdist(as.matrix(column_to_rownames(pe_wide3, "PlasmidID")), "jaccard", binary = T)

pe_clust <- hclust(plasmids_gene_dist, method = "average")
hclust_plasmids <- cutree(pe_clust, h = 0.5)
plasmid_clust <- data.table(plasmid = names(hclust_plasmids), PlasmidClust = hclust_plasmids)
potential_plasmids2 <- merge(potential_plasmids2, plasmid_clust, by.x = "PlasmidID", by.y = "plasmid", all.x = T)

bakta_im_plasmids <- merge(bakta_im_plasmids, potential_plasmids2[,list(genome, baktaContig, contig,  
                                                                        PlasmidClust)], by.x = c("Contig", "genome"), 
                           by.y = c("baktaContig", "genome"), all = T)

bakta_im_plasmids <- merge(bakta_im_plasmids, plasmid_eggnog[,list(PlasmidID, genome, `#query`, eggNOG_OGs, Description, 
                                                                   TopLevel_eggNOG, COG_category, COG_category2)], 
                           by.x = c("genome", "PlasmidID", "Locus Tag"), by.y = c("genome", "PlasmidID", "#query"), 
                           all.x = T)


cat_names <- data.table(COG_category = c(el_eggnog[nchar(COG_category)==1, sort(unique(COG_category))], "W", "Z"), 
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

bakta_im_plasmids <- merge(bakta_im_plasmids, cat_names, by = "COG_category", all.x = T)
bakta_im_plasmids[,COG_category3:=case_when(
  is.na(COG_category2)|COG_category2 == "-" ~ "None", 
  COG_category2=="Multiple/Other" ~ "Multiple/Other",
  TRUE ~ paste0(COG_category2, " ", CategoryName))]

cats_all <- bakta_im_plasmids[,sort(unique(COG_category3))]
colors_all <- c(brewer.pal(9, "Pastel1"), c(brewer.pal(3, "Pastel2")[2]))
names(colors_all) <- cats_all

## Make plasmid maps - Figures 5A, S5E-F
dir.create("example_plasmid_maps/")
source("plasmapR_modified.R")
make_save_plasmid_map <- function(plasmid_id, save_pdf = T, fillcolor = "COG", geneLab = "Product"){
  annots1 <- bakta_im_plasmids[PlasmidID==plasmid_id]
  if(geneLab == "Product"){
    setnames(annots1, "Product", "name")
  } else if(geneLab == "Description"){
    setnames(annots1, "Description", "name")
    annots1[,Description:=str_trunc(Description, width = 30)]
  }
  setnames(annots1, c("Start", "Stop"), c("start", "end"))
  #setnames(annots1, "Strand", "direction")
  annots1[,direction:=ifelse(Strand=="+", 1, -1)]
  #annots1[,type:="CDS"]
  if(fillcolor == "COG"){
    annots1[,type:=COG_category3]
  } else if(fillcolor == "Phylum"){
    annots1[,type:=PhylumLevel]
  } else if(fillcolor == "CorePrev"){
    annots1[,type:=CorePrev]
  } else {
    annots1[,type:=get(fillcolor)]
  }
  #annots1[is.na(type), type:="None"]
  annots1 <- unique(annots1[,list(name, type, start, end, direction)])
  annots1[,index:=2:(nrow(annots1)+1)]
  plot1 <- plot_plasmid2(annots1, bp = potential_plasmids2[PlasmidID==plasmid_id, len], name = paste0(
    gsub("contig", " contig", plasmid_id), "\nPlasmid Cluster ", potential_plasmids2[PlasmidID==plasmid_id, PlasmidClust]), label_wrap = 15)
  if(fillcolor == "COG"){
    plot1 <- plot1 +scale_fill_manual(values = colors_all, name = "COG Category", guide = "none")
  } else if(fillcolor == "Phylum"){
    plot1 <- plot1 +scale_fill_brewer(palette = "Set2", name = "eggNOG Phylum", guide = "none")
  }
  if(save_pdf){
    ggsave(plot1, file = paste0("example_plasmid_maps/", plasmid_id, ".pdf"), width = 9, height = 9)
  }
  return(plot1)
}

for(j in 1:length(potential_plasmids2[,PlasmidID])){
  make_save_plasmid_map(potential_plasmids2[j,PlasmidID])
}

## Write out potential plasmid seqs to file
for(i in 1:nrow(potential_plasmids2)){
  foo <- readDNAStringSet(paste0("elenta_fnas/", potential_plasmids2[i, genome], ".fna"))
  seq1 <- foo[names(foo) == potential_plasmids2[i, contig]]
  names(seq1) <- potential_plasmids2[i, PlasmidID]
  writeXStringSet(seq1, file = paste0("highConf_plasmids/", potential_plasmids2[i,PlasmidID], ".fna"))
}



## BLAST results from whole dataset - from blast_plasmids.sh
blast_results_fil <- list.files("plasmid_blast_results_rev/", full.names = T)
blast_results <- rbindlist(lapply(blast_results_fil, function(x){
  foo <- fread(x, header = F)
  return(foo)
}))
setnames(blast_results, paste0("V", 1:12), 
         c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
           "qend", "sstart", "send", "evalue", "bitscore"))
blast_results[,matchGenome:=gsub("__.*", "", sseqid)]
blast_results <- merge(blast_results, potential_plasmids2[,list(PlasmidID, len, PlasmidClust, baktaContig, contig)], by.x = "qseqid", by.y = "PlasmidID", all.x = T)

# Save top blast results that cover most of the sequence
blast_results[,matchRank:=length(bitscore)-rank(bitscore), by = list(qseqid, matchGenome)]
blast_results_sub <- blast_results[matchRank==0]
blast_results_sub[,totalLengthByContig:=sum(length), by = list(qseqid, sseqid)]
blast_results_sub[,clustRank:=length(bitscore)-rank(bitscore), by = list(sseqid, PlasmidClust)]
blast_results_sub[,FracLen:=length/len]
blast_results_sub <- blast_results_sub[totalLengthByContig > 1000|FracLen > 0.25]

blast_results_sub[,matchContig:=gsub(".*__", "", sseqid)]
blast_results_sub <- merge(blast_results_sub, contig_lengths, by.x = c("matchGenome", "matchContig"), by.y = c("genome", "contig"), all.x = T)
setnames(blast_results_sub, c("len.x", "len.y"), c("qseq_len", "sseq_len"))
blast_results_sub[,FracLenSSeq:=length/sseq_len]

blast_results_sub <- blast_results_sub[matchGenome != "GCF_900184265.1_PRJEB20940_genomic"] #timonensis

# remove matches that are shorter and lower qulaity than another hit for the same region
blast_results_sub[,internalMatch:=sapply(1:nrow(blast_results_sub), function(x){
  qi <- blast_results_sub[x, qseqid]
  si <- blast_results_sub[x, sseqid]
  ss <- blast_results_sub[x, sstart]
  se <- blast_results_sub[x, send]
  bits <- blast_results_sub[x, bitscore]
  foo <- blast_results_sub[sseqid==si ]
  if(nrow(foo) == 0) return(0) else {
    #Larger regions than the one in question
    if(ss < se){
      foo_sub <- foo[sstart <= ss & send >= se]
    } else { #Reverse complement
      foo_sub <- foo[sstart >= ss & send <= se]
    }
    foo_sub <- foo_sub[bitscore > bits] # This will also remove lower quality hits to the exact same region
    if(nrow(foo_sub) > 0) return(1) else return(0)
  }
})]
blast_results_sub <- blast_results_sub[internalMatch == 0]

## Now we have presence/absence of each cluster across all genomes
confident_present <- blast_results_sub[FracLen > 0.6 & length > 1000 & pident > 85]
confident_present_count <- confident_present[,list(length[which.max(length)], FracLen[which.max(length)], pident[which.max(length)]), by = list(matchGenome, matchContig, sseqid, PlasmidClust)]
confident_present_count <- merge(confident_present_count, el_metadata[,list(genome, Clust2, Continent, IsMAG)], by.x = "matchGenome", by.y = "genome")
common_clusts <- confident_present_count[,length(unique(matchGenome)), by=PlasmidClust][order(V1, decreasing = T)][V1 > 2, PlasmidClust]

genome_plasmid_count <- confident_present_count[PlasmidClust %in% common_clusts,length(unique(matchContig)), by = list(matchGenome, PlasmidClust)] %>% 
  dcast(matchGenome~PlasmidClust, value.var = "V1", fill = 0)
setnames(genome_plasmid_count, c("genome", paste0("PlasmidClust", names(genome_plasmid_count)[2:ncol(genome_plasmid_count)])))

# Add to metadata
el_metadata <- merge(el_metadata, genome_plasmid_count, by = "genome", all.x = T)

# How does prevalence in the whole dataset compare with prevalence in just the high-q? 
blast_count <- confident_present_count[,list(length(unique(matchContig)), median(V2), median(as.numeric(V1))), by = PlasmidClust]
setnames(blast_count, c("V1", "V2", "V3"), c("NumContigs", "MedianFracLen", "MedianBpBlast"))
blast_count[,Prevalence:=NumContigs/nrow(el_metadata)]
orig_count <- potential_plasmids2[,list(length(PlasmidID), median(as.numeric(len))), by = PlasmidClust]
setnames(orig_count, c("V1", "V2"), c("NumContigs", "MedianBpOrig"))

orig_genome_set <- el_metadata[genome %in% gsub("_circle.bed","", circle_map_files), genome]
orig_count[,Prevalence:=NumContigs/nrow(el_metadata)]
prev_compare <- merge(orig_count, blast_count, by = "PlasmidClust", all = T)
prev_compare[is.na(Prevalence.y), Prevalence.y:=0]

blast_count_isMAG <- confident_present_count[,list(length(unique(matchContig)), 
                                                   median(V2), median(as.numeric(V1))), 
                                             by = list(PlasmidClust, IsMAG)] %>% 
  dcast(PlasmidClust~IsMAG, value.var = c("V1", "V2", "V3"), fill = 0)
blast_count_isMAG[,PrevIsolate:=V1_FALSE/nrow(el_metadata[IsMAG==F])]
blast_count_isMAG[,PrevMAG:=V1_TRUE/nrow(el_metadata[IsMAG==T])]
blast_count_isMAG[,AvgSizeMAG:=ifelse(V3_TRUE==0, NA, V3_TRUE)]
blast_count_isMAG[,AvgSizeIsolate:=ifelse(V3_FALSE==0, NA, V3_FALSE)]

blast_count_isMAG[,labelConj:=case_when(
  PlasmidClust==2 ~"Conjugative", 
  PlasmidClust == 3 ~ "Cluster 2",
  PlasmidClust==7 ~ "Cluster 3",
  TRUE ~ "")]

# Figure S5A
isolate_mag_clust_compare <- ggplot(blast_count_isMAG, aes(x = PrevIsolate, y = PrevMAG)) + 
  geom_jitter(height = 0.005, width = 0.005, alpha = 0.8) + geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray") + 
  xlab("Prevalence in isolates") + ylab("Prevalence in MAGs") + xlim(0, 0.2) + ylim(0, 0.2) + 
  geom_text_repel(aes(label = labelConj), max.overlaps = 30)
ggsave(isolate_mag_clust_compare, file = "plasmids_isolateMagPrevCompare.pdf", width = 4, height = 3.8)

# Fig S5B-D
isolate_mag_size_compare <- ggplot(confident_present_count[PlasmidClust %in% c(2,3,7)], aes(x = ifelse(IsMAG, "MAGs", "Isolates"), y = V1)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.15, alpha = 0.7) + facet_wrap(~PlasmidClust, scales = "free") + xlab("") + ylab("Size of element (bp)") 

ggsave(isolate_mag_size_compare, file = "plasmids_clust237_isolateMAGSizeCompare.pdf", width = 6, height = 3.8)

el_metadata[,InOrigPlasmidSet:=ifelse(genome %in% orig_genome_set, 1, 0)]


pvals <- el_metadata %>% 
  pivot_longer(cols = starts_with("PlasmidClust"), names_to = "PlasmidClust", values_to = "Present") %>% 
  filter(PlasmidClust %in% paste0("PlasmidClust", common_clusts)) %>% 
  mutate(Present = ifelse(Present == 0 | is.na(Present), 0, 1)) %>% 
  group_by(PlasmidClust) %>% 
  summarize(OR = fisher.test(Clust2, Present)$estimate, PVal = fisher.test(Clust2, Present)$p.value) %>% 
  mutate(PAdj = p.adjust(PVal, method ="BH")) %>% 
  mutate(Signif = case_when(
    PAdj < 0.001 ~ "***", 
    PAdj < 0.01 ~ "**", 
    PAdj < 0.05 ~ "*", 
    TRUE ~ ""))

## Genome source check
el_metadata <- el_metadata %>% mutate(
  GenomeSource2 = case_when(
    GenomeSource.x %in% c("Harris Wang lab", "ImmunoMicrobiome") ~ "This study", 
    GenomeSource.x %in% c("NCBI", "Other") ~ "NCBI/Other",
    TRUE ~ GenomeSource.x
  )
)
pvals_source <- el_metadata %>% 
  pivot_longer(cols = starts_with("PlasmidClust"), names_to = "PlasmidClust", values_to = "Present") %>% 
  filter(PlasmidClust %in% paste0("PlasmidClust", c(2,3,7))) %>% 
  mutate(Present = ifelse(Present == 0 | is.na(Present), 0, 1)) %>% 
  group_by(PlasmidClust) %>% 
  summarize(OR = fisher.test(GenomeSource2, Present)$estimate, PVal = fisher.test(GenomeSource2, Present)$p.value) %>% 
  mutate(PAdj = p.adjust(PVal, method ="BH")) %>% 
  mutate(Signif = case_when(
    PAdj < 0.001 ~ "***", 
    PAdj < 0.01 ~ "**", 
    PAdj < 0.05 ~ "*", 
    TRUE ~ ""))

## Model
el_metadata <- el_metadata[,Continent2:=ifelse(Continent %in% c("Americas", "Asia", "Europe"), Continent, "Other/unknown")]

glm_res <- el_metadata %>% 
  pivot_longer(cols = starts_with("PlasmidClust"), names_to = "PlasmidClust", values_to = "Present") %>% 
  filter(PlasmidClust %in% paste0("PlasmidClust", common_clusts[1:3])) %>% 
  mutate(Present = ifelse(Present == 0 | is.na(Present), 0, 1)) %>% 
  group_by(PlasmidClust) %>% 
  reframe(broom::tidy(glm(Present~yearCompiled+factor(IsMAG)+Completeness+GenomeSource2+Continent2+Clust2, family = "binomial"))) %>% 
  group_by(term) %>% 
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  ungroup()
glm_res %>% filter(p.adj < 0.1)

glm_res <- el_metadata %>% 
  pivot_longer(cols = starts_with("PlasmidClust"), names_to = "PlasmidClust", values_to = "Present") %>% 
  filter(PlasmidClust %in% paste0("PlasmidClust", common_clusts[1:3])) %>% 
  mutate(Present = ifelse(Present == 0 | is.na(Present), 0, 1)) %>% 
  group_by(PlasmidClust) %>% 
  reframe(broom::tidy(glm(Present~factor(IsMAG)+Completeness+Continent2+Clust2, family = "binomial"))) %>% 
  group_by(term) %>% 
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  ungroup()
glm_res %>% filter(p.adj < 0.1)


## Get all genes, save for Supp Table
bakta_dat_sub <- data.table()
for(j in 1:nrow(confident_present)){
  foo <- bakta_dat[genome == confident_present[j, matchGenome] & Contig == confident_present[j, baktaContig.y] ]
  if(confident_present[j, sstart] < confident_present[j, send]){
    foo <- foo[Start >= confident_present[j, sstart] & Stop <= confident_present[j, send]]
  } else {
    foo <- foo[Start >= confident_present[j, send] & Stop <= confident_present[j, sstart]]
  }
  foo[,PlasmidClust:=confident_present[j, PlasmidClust]]
  bakta_dat_sub <- rbind(bakta_dat_sub, foo, fill = T)
}
bakta_dat_sub <- unique(bakta_dat_sub)
bakta_dat_sub <- merge(bakta_dat_sub, el_eggnog[,list(`#query`, seed_ortholog, eggNOG_OGs, TopLevel_eggNOG, 
                                                      COG_category, COG_category2, Description)], by.x = "Locus Tag", by.y = "#query", all.x = T)
bakta_dat_sub <- merge(bakta_dat_sub, po_plasmid[value != "*"], by.x = "Locus Tag", by.y = "value", all.x = T)

## Save info on Plasmid Clusters besides the conjugative
annots_save <- bakta_dat_sub[PlasmidClust %in% c(3,7)]
annots_save[,BaktaUniref50:=gsub("UniRef50_", "", str_extract(DbXrefs, pattern = "UniRef50_.[0-9A-Z]+"))]
annots_summ <- annots_save[,length(unique(genome)), by = list(BaktaUniref50, PlasmidClust)]
annots_genes <- unique(annots_save[,list(PlasmidClust, BaktaUniref50, Gene, Product, DbXrefs)])
annots_genes <- nest(annots_genes, Gene = Gene, Product = Product, DbXrefs = DbXrefs)
annots_genes[,GeneList:=sapply(Gene, function(x){
  foo <- unique(unlist(x)[unlist(x)!= ""])
  if(length(foo) == 0) return("") else return(paste0(foo, collapse = ","))
})]
annots_genes[,ProductList:=sapply(Product, function(x){
  foo <- unique(unlist(x)[unlist(x)!= ""])
  if(length(foo) == 0) return("") else return(paste0(foo, collapse = ","))
})]
annots_genes[,DbList:=sapply(DbXrefs, function(x){
  foo <- unique(unlist(x)[unlist(x)!= ""])
  if(length(foo) == 0) return("") else{
    foo <- unique(unlist(strsplit(foo, split = ", ")[[1]]))
    return(paste0(foo, collapse = ","))
  }
})]
annots_summ <- merge(annots_summ, annots_genes[,list(PlasmidClust, BaktaUniref50, GeneList, ProductList, DbList)], by = c("PlasmidClust", "BaktaUniref50"), all.x = T)

bakta_dat_sub[,length(unique(genome)), by = PlasmidClust][PlasmidClust %in% c(2,3,7)]
annots_summ[,ClustPrev:=ifelse(PlasmidClust==3, V1/35, V1/50)]
annots_summ[,ConservationLevel:=case_when(
  ClustPrev > 0.9 ~ "Nearly conserved (>90%)",
  TRUE ~ "<90% conserved"
)]

bakta_dat_sub[,BaktaUniref50:=gsub(", .*", "", gsub(".*UniRef:UniRef50_", "", DbXrefs))]
rep_genes <- bakta_dat_sub[(genome=="HS64-1A" & PlasmidClust==3)|(genome=="HS50-1A" & PlasmidClust==7), list(BaktaUniref50, `Locus Tag`, genome, PlasmidClust)]
annots_summ <- merge(annots_summ, rep_genes, by = c( "PlasmidClust", "BaktaUniref50"), all.x = T)
annots_summ[,LocusTag:=paste0(`Locus Tag`, collapse = ","), by = list(PlasmidClust, BaktaUniref50)]
annots_summ[,`Locus Tag`:=NULL]
annots_summ <- unique(annots_summ)
annots_summ[PlasmidClust == 3, PlasmidCluster:=2]
annots_summ[PlasmidClust == 7, PlasmidCluster:=3]

fwrite(annots_summ[,list(PlasmidCluster, genome, LocusTag, V1, ClustPrev, ConservationLevel, BaktaUniref50, GeneList, ProductList, DbList)] %>% arrange(PlasmidCluster, desc(V1)), 
       file = "plasmid_clust3_7_annotation_summary.csv")


###################### More detailed analysis of Conjugative plasmid genes

# First - reorient plasmids so they all have the same starting coordinate position
plasmid_seqs <- readDNAStringSet("highConf_plasmids/PlasmidCluster2/allPlasmidClust2.fna")

bakta_plasmid <- data.table()
for(j in 1:length(plasmid_seqs)){
  genome_id <- gsub("contig.*", "", names(plasmid_seqs)[j])
  contig <- gsub(".*contig_", "", names(plasmid_seqs)[j])
  foo <- fread(paste0("annots_out/bakta_elentaOnly/", genome_id, ".tsv"))
  plasmid_genes <- foo[`#Sequence Id`==paste0("contig_", contig)]
  plasmid_genes[,PlasmidID:=names(plasmid_seqs)[j]]
  bakta_plasmid <- rbind(bakta_plasmid, plasmid_genes)
}
names(plasmid_seqs)[!names(plasmid_seqs) %in% bakta_plasmid[grepl("Relaxase", Product), PlasmidID]]
names(plasmid_seqs)[!names(plasmid_seqs) %in% bakta_plasmid[grepl("Smc", Product), PlasmidID]]

starting_coordinate <- bakta_plasmid[Product=="DNA polymerase IV", list(PlasmidID, Start, Stop, Strand)]
starting_coordinate[,len:=nchar(plasmid_seqs)]
reoriented_seqs <- DNAStringSet()
for(j in 1:length(plasmid_seqs)){
  if(starting_coordinate[j, Strand]=="+"){
    new_start <- starting_coordinate[j, Start]
    new_seq <- DNAString(substr(plasmid_seqs[j], new_start, nchar(plasmid_seqs[j])))
    new_seq <- append(new_seq, DNAString(substr(plasmid_seqs[j], 1, new_start-1)))
  } else {
    new_start <- starting_coordinate[j, Stop]
    new_seq <- reverseComplement(DNAString(substr(plasmid_seqs[j], 1, new_start)))
    new_seq <- append(new_seq, reverseComplement(DNAString(substr(plasmid_seqs[j], new_start+1, nchar(plasmid_seqs[j])))))
  }
  reoriented_seqs <- append(reoriented_seqs, DNAStringSet(new_seq))
}
names(reoriented_seqs) <- names(plasmid_seqs)

aligned_seqs <- AlignSeqs(reoriented_seqs)
nchar(aligned_seqs)[6]-50980
## Fix HS21-3A - move 93bp from end to beginning
substr(aligned_seqs[6], 50981, nchar(aligned_seqs)[6])
new21 <- append(DNAString(substr(aligned_seqs[6],  50981, nchar(aligned_seqs)[6])), DNAString(substr(reoriented_seqs[6], 1, nchar(reoriented_seqs[6])-93)))
reoriented_seqs[[6]] <- new21
aligned_seqs <- AlignSeqs(reoriented_seqs)
#BrowseSeqs(aligned_seqs)

writeXStringSet(reoriented_seqs, file = "highConf_plasmids/PlasmidCluster2/PlasmidCluster2_reoriented.fna")

## Update gene coordinates based on reoriented sequences
reoriented_coords <- data.table()
for(j in 1:length(plasmid_seqs)){
  foo <- bakta_plasmid[PlasmidID==names(plasmid_seqs)[j]]
  if(starting_coordinate[j, Strand]=="+"){
    offset <- starting_coordinate[j, Start]
    foo[,NewStart:=Start+1-offset]
    foo[,NewStop:=Stop+1-offset]
    foo[NewStart < 0, NewStart:=starting_coordinate[j, len]+NewStart]
    foo[NewStop < 0, NewStop:=starting_coordinate[j, len]+NewStop]
    foo[,NewStrand:=Strand]
  } else {
    # First just reverse the coordinates, then move backwards
    foo[,NewStrand:=ifelse(Strand=="-", "+", "-")]
    foo[,NewStart1:=starting_coordinate[j, len]+1-Stop]
    foo[,NewStop1:=starting_coordinate[j, len]+1-Start]
    offset <-foo[Product=="DNA polymerase IV", NewStart1]
    foo[,NewStart:=NewStart1+1-offset]
    foo[,NewStop:=NewStop1+1-offset]
    foo[NewStart < 0, NewStart:=starting_coordinate[j, len]+NewStart]
    foo[NewStop < 0, NewStop:=starting_coordinate[j, len]+NewStop]
  }
  reoriented_coords <- rbind(reoriented_coords, 
                             foo[,list(PlasmidID, Start, Stop, Strand, NewStart, NewStop, NewStrand, `Locus Tag`, Gene, Product)], fill = T)
}
setnames(reoriented_coords, c("Start", "Stop", "Strand", "NewStart", "NewStop", "NewStrand"), 
         c("OldStart", "OldStop", "OldStrand", "Start", "Stop", "Strand"))
#Fix coordinates for HS21
reoriented_coords[PlasmidID=="HS21-3Acontig_2", Start:=Start+93]
reoriented_coords[PlasmidID=="HS21-3Acontig_2", Stop:=Stop+93]
reoriented_coords[PlasmidID=="HS21-3Acontig_2" & Product == "DNA polymerase IV", Start:=1]

### Updated plotting
reoriented_coords <- merge(reoriented_coords, bakta_im_plasmids[,list(`Locus Tag`, COG_category3)], by = "Locus Tag")
for(j in 1:length(starting_coordinate[,PlasmidID])){
  make_save_plasmid_map(starting_coordinate[j,PlasmidID], gene_tab = reoriented_coords[PlasmidID == starting_coordinate[j, PlasmidID]])
}

# Get new annotations, extract adhesin genes (COG4932) (+compare with old annots)
new_annots <- fread("plasmid_cluster2/bakta_out/PlasmidsClust2.tsv")
id_maps <- fread(cmd = "grep orig-id plasmid_cluster2/bakta_out/PlasmidsClust2.log", header = F, sep = " ")
id_maps[,PlasmidID:=gsub("orig-id=", "", gsub(",$", "", V11))]
id_maps[,baktaName:=gsub("id=", "", gsub(",$", "", V10))]
new_annots <- merge(new_annots, id_maps[,list(PlasmidID, baktaName)], by.x = "#Sequence Id", by.y = "baktaName")

plasmid_seqs <- readDNAStringSet("plasmid_cluster2/PlasmidCluster2_reoriented.fna")

## Do with blast results
# Get plasmids that were identified based on BLAST and add these to the gene annotations
clust2_contigs <- confident_present[PlasmidClust==2]

all_clust2_seqs <- DNAStringSet()
contigs_set <- clust2_contigs[,sort(unique(sseqid))]
for(j in 1:length(contigs_set)){
  geno <- gsub("__.*", "", contigs_set[j])
  cont <- gsub(".*__", "", contigs_set[j])
  if(!geno %in% gsub("contig.*", "", names(plasmid_seqs))){ # Only ones that we don't already have
    foo <- readDNAStringSet(paste0("eggerthellaceae_collection/fnas/", geno, ".fna"))
    foo <- foo[gsub(" .*", "", names(foo)) == cont]
    names(foo) <- paste0(geno, "__", clust2_contigs[matchContig==cont, unique(baktaContig.y)]) # Switch to bakta contigs
    all_clust2_seqs <- append(all_clust2_seqs, foo)
  }
}

bakta_plasmid_add <- data.table()
for(j in 1:length(all_clust2_seqs)){
  genome_id <- gsub("__contig.*", "", names(all_clust2_seqs)[j])
  contig <- gsub(".*contig_", "", names(all_clust2_seqs)[j])
  plasmid_genes <- bakta_dat[genome == genome_id & Contig==paste0("contig_", contig)]
  plasmid_genes[,PlasmidID:=names(all_clust2_seqs)[j]]
  bakta_plasmid_add <- rbind(bakta_plasmid_add, plasmid_genes)
}
starting_coordinate <- bakta_plasmid_add[Product=="DNA polymerase IV", list(PlasmidID, Start, Stop, Strand)]
starting_coordinate[,len:=nchar(all_clust2_seqs)]
reoriented_seqs <- DNAStringSet()
for(j in 1:length(all_clust2_seqs)){
  if(starting_coordinate[j, Strand]=="+"){
    new_start <- starting_coordinate[j, Start]
    new_seq <- DNAString(substr(all_clust2_seqs[j], new_start, nchar(all_clust2_seqs[j])))
    new_seq <- append(new_seq, DNAString(substr(all_clust2_seqs[j], 1, new_start-1)))
  } else {
    new_start <- starting_coordinate[j, Stop]
    new_seq <- reverseComplement(DNAString(substr(all_clust2_seqs[j], 1, new_start)))
    new_seq <- append(new_seq, reverseComplement(DNAString(substr(all_clust2_seqs[j], new_start+1, nchar(all_clust2_seqs[j])))))
  }
  reoriented_seqs <- append(reoriented_seqs, DNAStringSet(new_seq))
}
names(reoriented_seqs) <- names(all_clust2_seqs)


combined_clust2 <- append(plasmid_seqs, reoriented_seqs)
aligned_seqs <- AlignSeqs(combined_clust2)

## Save new compilation of all sequences with reoriented coordinates
writeXStringSet(combined_clust2, file = "plasmid_cluster2/AllPlasmidCluster2_reoriented.fna")

## Fix annotation coordinates also
reoriented_coords <- data.table()
for(j in 1:length(reoriented_seqs)){
  foo <- bakta_plasmid_add[PlasmidID==names(reoriented_seqs)[j]]
  if(starting_coordinate[j, Strand]=="+"){
    offset <- starting_coordinate[j, Start]
    foo[,NewStart:=Start+1-offset]
    foo[,NewStop:=Stop+1-offset]
    foo[NewStart < 0, NewStart:=starting_coordinate[j, len]+NewStart]
    foo[NewStop < 0, NewStop:=starting_coordinate[j, len]+NewStop]
    foo[,NewStrand:=Strand]
  } else {
    # First just reverse the coordinates, then move backwards
    foo[,NewStrand:=ifelse(Strand=="-", "+", "-")]
    foo[,NewStart1:=starting_coordinate[j, len]+1-Stop]
    foo[,NewStop1:=starting_coordinate[j, len]+1-Start]
    offset <-foo[Product=="DNA polymerase IV", NewStart1]
    foo[,NewStart:=NewStart1+1-offset]
    foo[,NewStop:=NewStop1+1-offset]
    foo[NewStart < 0, NewStart:=starting_coordinate[j, len]+NewStart]
    foo[NewStop < 0, NewStop:=starting_coordinate[j, len]+NewStop]
  }
  reoriented_coords <- rbind(reoriented_coords, 
                             foo[,list(PlasmidID, Start, Stop, Strand, NewStart, NewStop, NewStrand, `Locus Tag`, Gene, Product)], fill = T)
}
setnames(reoriented_coords, c("Start", "Stop", "Strand", "NewStart", "NewStop", "NewStrand"), 
         c("OldStart", "OldStop", "OldStrand", "Start", "Stop", "Strand"))
reoriented_coords[,genome:=gsub("__.*", "", PlasmidID)]
reoriented_coords <- merge(reoriented_coords, el_eggnog[,list(genome, `#query`, eggNOG_OGs, Description, 
                                                              TopLevel_eggNOG, COG_category, COG_category2)], 
                           by.x = c("genome",  "Locus Tag"), by.y = c("genome","#query"), 
                           all.x = T)

## Add COG categories
cat_names <- data.table(COG_category = c(el_eggnog[nchar(COG_category)==1, sort(unique(COG_category))], "W", "Z"), 
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

reoriented_coords <- merge(reoriented_coords, cat_names, by = "COG_category", all.x = T)
reoriented_coords[,COG_category3:=case_when(
  is.na(COG_category2)|COG_category2 == "-" ~ "None", 
  COG_category2=="Multiple/Other" ~ "Multiple/Other",
  TRUE ~ paste0(COG_category2, " ", CategoryName))]


### Use finalized dataset
new_annots_all <- fread("plasmid_cluster2/bakta_out/AllPlasmidsClust2.tsv")
plasmid_seqs <- readDNAStringSet("plasmid_cluster2/AllPlasmidCluster2_reoriented.fna")

id_maps <- fread(cmd = "grep orig-id plasmid_cluster2/bakta_out/AllPlasmidsClust2.log", header = F, sep = " ")
id_maps[,PlasmidID:=gsub("orig-id=", "", gsub(",$", "", V11))]
id_maps[,baktaName:=gsub("id=", "", gsub(",$", "", V10))]
new_annots_all <- merge(new_annots_all, id_maps[,list(PlasmidID, baktaName)], by.x = "#Sequence Id", by.y = "baktaName")
plasmids_all <- new_annots_all[,sort(unique(PlasmidID))]

## Review variation in gene annotations between the plasmids
new_annots_all[,BaktaUniref50:=gsub(", .*", "", gsub(".*UniRef:UniRef50_", "", DbXrefs))]
new_annots_all[,length(unique(PlasmidID)), by = list(Product, BaktaUniref50)][order(V1, decreasing = T)]
new_annots_all[,length(unique(PlasmidID)), by = list(BaktaUniref50)][order(V1, decreasing = T)][V1 > 33]
fully_present <- new_annots_all[,length(unique(PlasmidID)), by = list(BaktaUniref50)][V1==66, BaktaUniref50]
mostly_present <- new_annots_all[,length(unique(PlasmidID)), by = list(BaktaUniref50)][order(V1, decreasing = T)][V1 > 59 & V1 < 66, BaktaUniref50]
new_annots_all[BaktaUniref50 %in% fully_present, unique(Product)]
new_annots_all[BaktaUniref50 %in% mostly_present, unique(Product)]
new_annots_all[,ConservationLevel:=factor(case_when(
  BaktaUniref50 %in% fully_present ~ "Always present",
  BaktaUniref50 %in% mostly_present ~ "Nearly conserved (>90%)",
  TRUE ~ "<90% conserved"
), levels = c("Always present", "Nearly conserved (>90%)", "<90% conserved"))]

new_annots_all[grepl("PrgI", Product),length(unique(PlasmidID)), by = list(BaktaUniref50)]
plasmids_all[sapply(plasmids_all, function(x){
  nrow(new_annots_all[PlasmidID==x & BaktaUniref50 %in% mostly_present]) < 18
})]
# Mainly 3 instances missing conjugation proteins, some variation in tox-antitox, helicase? 
new_annots_all[BaktaUniref50 %in% mostly_present[!mostly_present %in% new_annots_all[PlasmidID=="HS23-3A__contig_4", BaktaUniref50]], unique(Product)]
new_annots_all[BaktaUniref50 %in% mostly_present[!mostly_present %in% new_annots_all[PlasmidID=="HS108-1Acontig_5", BaktaUniref50]], unique(Product)]
new_annots_all[BaktaUniref50 %in% mostly_present[!mostly_present %in% new_annots_all[PlasmidID=="GCA_958412325.1_ERR10960891_bin.8_MetaWRAP_v1.3_MAG_genomic__contig_24", BaktaUniref50]], unique(Product)]

# Figure 5B - final version of prevalence graph
new_annots_all[,genome:=gsub("__$", "", gsub("contig.*", "", PlasmidID))]
el_metadata[,HasPlasmidClust2:=ifelse(genome %in% new_annots_all[,genome], 1, 0)]
el_metadata[,Continent2:=factor(ifelse(Continent %in% c("Americas", "Asia", "Europe"), Continent, "Other/unknown"), 
                                levels = c("Other/unknown", "Europe", "Asia", "Americas"))]
el_metadata <- el_metadata %>% mutate(
  GenomeSource2 = case_when(
    GenomeSource.x %in% c("Harris Wang lab", "ImmunoMicrobiome") ~ "This study", 
    GenomeSource.x %in% c("NCBI", "Other") ~ "NCBI/Other",
    TRUE ~ GenomeSource.x
  )
)

genome_list <- gsub("__.*$", "", names(plasmid_seqs))
genome_list <- gsub("contig.*$", "", genome_list)
el_metadata[,HasPlasmidClust2:=ifelse(genome %in% genome_list, 1, 0)]
clust_prev <- ggplot(el_metadata, aes(x = ifelse(Clust2==2, "B", "A"), fill = factor(HasPlasmidClust2))) + geom_bar(color = "black", position = "fill") + coord_flip() + 
  facet_wrap(~ifelse(IsMAG, "MAGs", "Isolates")) + theme_classic2() + xlab("Genome group") + ylab("Fraction of genomes with plasmid detected") +
  scale_fill_manual(values = c("white", "darkblue"), guide = "none")
ggsave(clust_prev, file = "conjplasmid_prev_clust.pdf", width = 5, height = 2.7)

el_metadata[,fisher.test(Clust2, HasPlasmidClust2)]
el_metadata[,summary(glm(HasPlasmidClust2~IsMAG+Clust2, family = "binomial"))] ## Not signif

el_metadata[,summary(glm(HasPlasmidClust2~IsMAG+Continent2, family = "binomial"))] # Asia is borderline signif in this formulation

## Save general annotation info for supp table S5

## Plasmid map diagram with extra annotations - Figure 5A
plasmid_manual <- read_xlsx("ConjugativePlasmid_geneAnnotations.xlsx", sheet = 3) %>% data.table()

new_annots_all_sub <- new_annots_all[grepl("106-3A", PlasmidID)]
new_annots_all_sub <- merge(new_annots_all_sub, plasmid_manual[,list(`Locus Tag`, Category)], by = "Locus Tag" ,all.x = T)
new_annots_all_sub[is.na(Category), Category:=""]
new_annots_all_sub[,arrowColor:=ConservationLevel]
plasmid_map_focus <- make_save_plasmid_map("HS106-3Acontig_3", fillcolor = "Category") + 
  scale_fill_manual(values = c("gray93", brewer.pal(8, "Pastel2")), 
                    guide = guide_legend(override.aes = list(label = "")), name = "") +
  scale_color_manual(values = c("black", "gray35", "gray70"), 
                     guide = guide_legend(override.aes = list(fill = "white")), name = "Gene conservation")#+ 
#  theme(legend.key.size = unit(0.2, "in"))
ggsave(plasmid_map_focus, file = "plasmid2ManualAnnots_HS106-3a.pdf", width = 8, height = 8)

## Check Conservation - 100% or not
plasmid2_annots_summary <- new_annots_all[,length(unique(genome)), by = list(BaktaUniref50, ConservationLevel)]
plasmid2_genes <- unique(new_annots_all[,list(BaktaUniref50, Gene, Product, DbXrefs)])
plasmid2_genes <- nest(plasmid2_genes, Gene = Gene, Product = Product, DbXrefs = DbXrefs)
plasmid2_genes[,GeneList:=sapply(Gene, function(x){
  foo <- unique(unlist(x)[unlist(x)!= ""])
  if(length(foo) == 0) return("") else return(paste0(foo, collapse = ","))
})]
plasmid2_genes[,ProductList:=sapply(Product, function(x){
  foo <- unique(unlist(x)[unlist(x)!= ""])
  if(length(foo) == 0) return("") else return(paste0(foo, collapse = ","))
})]
plasmid2_genes[,DbList:=sapply(DbXrefs, function(x){
  foo <- unique(unlist(x)[unlist(x)!= ""])
  if(length(foo) == 0) return("") else{
    foo <- unique(unlist(strsplit(foo, split = ", ")[[1]]))
    return(paste0(foo, collapse = ","))
  }
})]
plasmid2_annots_summary <- merge(plasmid2_annots_summary, plasmid2_genes[,list(BaktaUniref50, GeneList, ProductList, DbList)], by = "BaktaUniref50", all.x = T)

## Combine with original annots - save for Table S5
bakta_plasmid[,genome:=gsub("contig.*", "", PlasmidID)]
bakta_plasmid[,BaktaUniref50:=gsub(", .*", "", gsub(".*UniRef:UniRef50_", "", DbXrefs))]
rep_genes <- bakta_plasmid[genome=="HS106-3A", list(BaktaUniref50, `Locus Tag`)]
plasmid2_annots_summary <- merge(plasmid2_annots_summary, rep_genes, by = "BaktaUniref50", all.x = T)
plasmid2_annots_summary[,LocusTag:=paste0(`Locus Tag`, collapse = ","), by = BaktaUniref50]
plasmid2_annots_summary[,`Locus Tag`:=NULL]
plasmid2_annots_summary <- unique(plasmid2_annots_summary)
plasmid2_annots_summary[LocusTag != "NA", genome:="HS106-3A"]

plasmid2_annots_summary[,PlasmidCluster:=1] #"Plasmid Cluster 1 (Conjugative)"]
fwrite(plasmid2_annots_summary[,list(PlasmidCluster, LocusTag, V1, ConservationLevel,BaktaUniref50, GeneList, ProductList, DbList)][order(V1, decreasing = T)], 
       file = "plasmid_clust2_annotation_summary.csv")


## Adhesin gene analysis with full dataset - Figure S5
adhesin_genes <- DNAStringSet()
for(j in 1:length(plasmids_all)){
  foo <- new_annots_all[PlasmidID==names(plasmid_seqs)[j] & (grepl("Prealbumin", Product)|grepl("SpaA", Product))]
  if(nrow(foo) ==1){
    geneseq <- DNAStringSet(substr(plasmid_seqs[j], foo[,Start], foo[,Stop]))
    names(geneseq) <- paste0(names(plasmid_seqs)[j], "__", foo[,Product])
    adhesin_genes <- append(adhesin_genes, geneseq)
  } else if(nrow(foo) > 1){
    for(k in 1:nrow(foo)){
      geneseq <- DNAStringSet(substr(plasmid_seqs[j], foo[k,Start], foo[k,Stop]))
      names(geneseq) <- paste0(names(plasmid_seqs)[j], "__", foo[k,Product], "__", k)
      adhesin_genes <- append(adhesin_genes, geneseq)
    }
  }
}
adhesin_prot <- translate(reverseComplement(adhesin_genes))
writeXStringSet(adhesin_prot, file = "plasmid_cluster2/allAdhesinSeqs.faa")
#Remove stop codons from file
system("sed -i -e 's/\\*//g' plasmid_cluster2/allAdhesinSeqs.faa")

aligned_prots <- AlignSeqs(adhesin_prot)

consenseq <- ConsensusSequence(aligned_prots, threshold = 0.3)

## Interpro domain annotations for these seqs
interpro_results <- fread("plasmid_cluster2/interpro_out/allAdhesinSeqs.faa.tsv")
interpro_results[,V15:=NULL]
pfam_tmm_results <- interpro_results[V4 %in% c("Pfam", "TMHMM")]
pfam_tmm_results[,genome:=gsub("__$", "", gsub("contig.*", "", V1))]
pfam_tmm_results <- merge(pfam_tmm_results, el_metadata[,list(genome, Continent2, Country, GenomeSource2, Max_Contig_Length, Contig_N50, n_contigs, gender, age_category, IsMAG)], by = "genome")

###### Figure S5G - Prevalence of Adhesin genes across both the plasmid and other instances 
adhesin_cog <- "COG4932"
adh_eg <- el_eggnog[grepl(adhesin_cog, eggNOG_OGs)]
adh_eg <- merge(adh_eg, bakta_dat, by.x = c("genome", "#query"), by.y = c("genome", "Locus Tag"))
adh_eg[,length(unique(`#query`)), by = genome][,table(V1)] # up to 9 per genome
#ComplexHeatmap::Heatmap(column_to_rownames(dcast(mutate(adh_eg, Present = 1), genome~eggNOG_OGs, value.var = "Present", fill = 0, fun.aggregate = sum), "genome")) #ggplot(adh_eg, aes(x = genome, y = eggNOG_OGs, fill = 1)) + geom_tile()

adh_eg[,length(unique(eggNOG_OGs))] # 22
interpro_adhesins <- rbindlist(lapply(el_metadata[,genome], function(x){
  foo <- fread(paste0("interpro_out_v2tmhmm/", x, "/", x, ".tsv"), quote = "", header = F)
  foo[,V15:=NULL]
  foo <- foo[V1 %in% adh_eg[,`#query`]]
  foo[,genome:=x]
  return(foo)
}))
pfam_adhesins_all <- interpro_adhesins[V4 %in% c("Pfam", "TMHMM")]
pfam_adhesins_all <- merge(pfam_adhesins_all, adh_eg[,list(`#query`, eggNOG_OGs, Description)], by.x = "V1", by.y = "#query", all.x = T)

rep_genes <- pfam_adhesins_all[genome %in% el_metadata[IsMAG==F, genome],list(V1[which.max(V3)], max(V3)), by = eggNOG_OGs]

num_genomes <- adh_eg[,length(unique(genome)), by = eggNOG_OGs][order(V1, decreasing = T)]
num_genomes <- num_genomes[eggNOG_OGs %in% rep_genes[,eggNOG_OGs]]
#genes have to appear at least twice
rep_genes <- rep_genes[eggNOG_OGs %in% num_genomes[V1 > 1, eggNOG_OGs]]
num_genomes <- num_genomes[eggNOG_OGs %in% rep_genes[,eggNOG_OGs]]
adh_domain_plot <- ggplot(pfam_adhesins_all[V1 %in% rep_genes[,V1]]) + 
  geom_segment(aes(x= 1, xend = V3, y = factor(eggNOG_OGs, levels =rev(num_genomes[,eggNOG_OGs]))), color = "black") +
  geom_segment(aes(x = V7, xend = V8, color = paste0(V6, " (", V5, ")"), y = factor(eggNOG_OGs, levels =rev(num_genomes[,eggNOG_OGs]))), linewidth = 2) + 
  #facet_wrap(~genome, scales = "free_y") +
  scale_color_manual(values = c(brewer.pal(8, "Set1"), brewer.pal(6, "Set2")), name = "Protein domain", 
                     guide = guide_legend(ncol = 2)) + 
  theme(legend.position = "bottom", axis.text.y = element_blank(), strip.text = element_text(size = 3), axis.title.y = element_blank()) +
  xlab("Amino acid coordinate") 

num_genomes_plot <- ggplot(num_genomes, aes(y = factor(eggNOG_OGs, levels = rev(num_genomes[,eggNOG_OGs])), 
                                            x = V1)) + 
  geom_bar(stat = "identity") +  scale_x_reverse()+
  scale_y_discrete(labels = gsub(".*,", "", rev(num_genomes[,eggNOG_OGs])), name = "Adhesin gene family", 
                   position = "right") +
  xlab("Number of genomes") 

annotate_tab <- data.table(eggNOG_OGs = c("COG4932@1|root,COG4932@2|Bacteria,2GKV5@201174|Actinobacteria,4CUR8@84998|Coriobacteriia", 
                                          "COG4932@1|root,COG4932@2|Bacteria,2HSRZ@201174|Actinobacteria,4CVDC@84998|Coriobacteriia",
                                          "COG4932@1|root,COG4932@2|Bacteria,1TS9Q@1239|Firmicutes,24AX1@186801|Clostridia,27WEU@189330|Dorea"), 
                           Label = c("Assoc with colonization in Bisanz 2020", "On Conj Plasmid", "On Conj Plasmid"))
final_adh_plot <- num_genomes_plot+adh_domain_plot+
  geom_label(data = annotate_tab, x = 2200, size = 3, aes(label = Label, 
                                                          y = factor(eggNOG_OGs, levels =rev(num_genomes[,eggNOG_OGs])))) + plot_layout(widths = c(1,2)) 
ggsave(final_adh_plot, file = "adhesin_variation_plot.pdf", width = 8, height = 6)

#####  Figure 5C - Analyze Plasmids across strains from same host
plasmid_seqs <- readDNAStringSet("plasmid_cluster2/AllPlasmidCluster2_reoriented.fna")
full_align <- AlignSeqs(plasmid_seqs)

im_seqs <- plasmid_seqs[grepl("HS", names(plasmid_seqs))]
im_hosts <- gsub("-.*", "", names(im_seqs))
host_dups <- names(table(im_hosts)[table(im_hosts) > 1])  # HS33, HS55, HS75
shared_sub <- im_seqs[sapply(names(im_seqs), function(y){ 
  any(str_detect(pattern = host_dups, y))
})]

#Calculate alignment distances between all plasmids
dists_all_plasmids <- DistanceMatrix(full_align, includeTerminalGaps = T, correction = "JC69")
dist_mat_compare <- dists_all_plasmids %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("Plasmid1") %>% 
  data.table() %>% melt(id.var = "Plasmid1", variable.name = "Plasmid2", variable.factor = F) %>% 
  filter(Plasmid1 < Plasmid2)
dist_mat_compare[grepl("HS", Plasmid1),Host1:=gsub("-.*", "", Plasmid1)]
dist_mat_compare[grepl("HS", Plasmid2),Host2:=gsub("-.*", "", Plasmid2)]
dist_mat_compare[!is.na(Host1) & !is.na(Host2),SameHost:=ifelse(Host1 == Host2, Host1, NA)]

## Calculate ANI between plasmids
ani_noIndels <- function(seqID1, seqID2, seqs_all){
  seqpair <- seqs_all[names(seqs_all) %in% c(seqID1, seqID2)]
  aligned_sub1 <- AlignSeqs(seqpair)
  consenmat1 <- consensusMatrix(aligned_sub1)
  indels <- which(apply(consenmat1, 2, function(x){
    x[names(x)=="-"] == 1
  }))
  print(indels)
  if(length(indels) > 0){
    consenmat_noIndels <- consenmat1[,-indels]
  } else {
    consenmat_noIndels <- consenmat1
  }
  print(paste0("Length of alignmnet, indels removed: ", dim(consenmat_noIndels)[2]))
  diffs1 <- which(apply(consenmat_noIndels, 2, function(x){
    any(!x %in% c(0, 2))
  }))
  ani_comp <- 1-length(diffs1)/ncol(consenmat_noIndels)
  return(ani_comp)
}
pairwise_anis <- data.table(expand_grid(plasmids_all, plasmids_all))
setnames(pairwise_anis, c("P1", "P2"))
pairwise_anis <- pairwise_anis[P1 < P2]
pairwise_anis[,ANI:=sapply(1:nrow(pairwise_anis), function(y){
  ani_noIndels(pairwise_anis[y, P1], pairwise_anis[y, P2], plasmid_seqs)
})]  
pairwise_anis[,hist(ANI, breaks = 200, xlim = c(0.97, 1.01))]
pairwise_anis[P1 %in% names(shared_sub) & P2 %in% names(shared_sub)]
fwrite(pairwise_anis, file = "plasmid_cluster2/pairwise_anis_all_plasmids.tsv", sep = "\t")
pairwise_anis <- fread("plasmid_cluster2/pairwise_anis_all_plasmids.tsv")

nrow(pairwise_anis[ANI >= 0.9998248])/nrow(pairwise_anis) #HS33
nrow(pairwise_anis[ANI >= 0.9997153])/nrow(pairwise_anis) #HS55
nrow(pairwise_anis[ANI ==1])/nrow(pairwise_anis) #HS75
nrow(pairwise_anis[ANI ==1])
pairwise_anis[ANI==1]
plasmids_all[!plasmids_all %in% pairwise_anis[,P1]]

pairwise_anis[((P1 == "HS75-3Acontig_4")|(P2 == "HS75-3Acontig_4")) & ANI != 1][order(ANI, decreasing = T)]
pairwise_anis[((P1 == "HS75-3Bcontig_2")|(P2 == "HS75-3Bcontig_2")) & ANI != 1][order(ANI, decreasing = T)]

ani_noIndels("HS75-3Bcontig_2", "GCF_015558765.1_ASM1555876v1_genomic__contig_12", plasmid_seqs)
45674*0.9999343 # only differ by 3 positions

align75 <- AlignSeqs(plasmid_seqs[c("HS75-3Acontig_4", "HS75-3Bcontig_2")])
consen75 <- consensusMatrix(align75)
which(apply(consen75, 2, max) != 2)
consen75[,33913] # Single insertion in an intergenic region

full_anis <- fread("pyani_elenta/ANIm_percentage_identity.tab") %>% 
  melt(id.var = "V1", variable.name ="V2", variable.factor = F) %>% 
  filter(V1 < V2)
full_anis[grepl("HS33", V1) & grepl("HS33", V2)] # 2A and 3A.- 98.7%
full_anis[grepl("HS55", V1) & grepl("HS55", V2)] # 2A and 3A.- 98.8%
full_anis[grepl("HS75", V1) & grepl("HS75", V2)] # 3A and 3B.- 98.4%

dist_mat_compare[,genome1:=gsub("__$", "", gsub("contig.*", "", Plasmid1))]
dist_mat_compare[,genome2:=gsub("__$", "", gsub("contig.*", "", Plasmid2))]
dist_mat_compare <- merge(dist_mat_compare, full_anis, by.x = c("genome1", "genome2"), by.y = c("V1", "V2"))
dist_mat_compare[,ANIdist:=1-value.y]

# Figure 5C
plasmid_samehost_dat <- ggplot(dist_mat_compare, aes(x = ANIdist, y = value.x)) + geom_point(alpha = 0.4, color = "darkgray") + 
  geom_point(data = dist_mat_compare[!is.na(SameHost)], aes(fill = SameHost), shape = 21, color = "black", size = 2.2) +
  scale_fill_brewer(palette = "Dark2", guide = "none") + xlab("Whole-genome ANI divergence") + 
  ylab("Plasmid sequence divergence (Jukes-Cantor)") + 
  geom_label_repel(data = dist_mat_compare[!is.na(SameHost)], aes(label = paste0(gsub("contig", "", Plasmid1), "/", gsub("contig", "", Plasmid2))), size = 2.2)
ggsave(plasmid_samehost_dat, file = "final_genome_set_2024-07/plasmid2_compare_plasmidDist_aniDist.pdf", width = 3.8, height = 3.8)
dist_mat_compare[,cor.test(ANIdist, value.x, method = "spearman")] ## Not significant

## Remove gaps present in greater than 1/2 of sequences
consen_mat <- consensusMatrix(full_align)
mostly_gaps <- which(colSums(consen_mat[5:18,]) >= 33)
align_mat <- as.matrix(full_align)
align_mat_sub <- align_mat[,-mostly_gaps]
full_align_sub <- DNAStringSet()
for(j in 1:nrow(align_mat_sub)){
  string1 <- DNAStringSet(paste0(align_mat_sub[j,], collapse = ""))
  full_align_sub <- append(full_align_sub, string1)
}
names(full_align_sub) <- row.names(align_mat_sub) 
full_plasmid_phylo <- Treeline(full_align_sub)

WriteDendrogram(full_plasmid_phylo, file = "conjplasmid_phylogeny.tre")

# Convert to phylo type
plasmid_phy <- read.tree("conjplasmid_phylogeny.tre")
plasmid_phy <- unroot(plasmid_phy)

# Bootstrap func that can do unrooted- use boot.phylo with a function where it finds the Treeline tree, writes and reads to temp file to get it in the right class, then compiles
get_tree_func <- function(x){
  string_set <- DNAStringSet()
  for(j in 1:nrow(x)){
    string1 <- DNAStringSet(paste0(x[j,], collapse = ""))
    string_set <- append(string_set, string1)
  }
  names(string_set) <- row.names(x) # Make sure the names are always the same
  tree1 <- Treeline(string_set)
  WriteDendrogram(tree1, file = "temp_dendro.tre")
  tree1p <- read.tree("temp_dendro.tre")
  tree1p <- unroot(tree1p)
  return(tree1p)
}

test_boot <- boot.phylo(plasmid_phy, align_mat_sub, FUN = get_tree_func, B = 100, trees = T, rooted = F)
plasmid_phy$node.label <- test_boot$BP

plasmid_phy$tip.label <- gsub("'", "", plasmid_phy$tip.label, fixed = T)
plasmid_tab <- data.table(Plasmid = plasmid_phy$tip.label, genome = gsub("__$", "", gsub("contig.*", "", plasmid_phy$tip.label)))
plasmid_tab <- merge(plasmid_tab, el_metadata[,list(genome, Continent, Clust2)], by = "genome")
plasmid_tab[,lab1:=ifelse(grepl("HS55", Plasmid)|grepl("HS75", Plasmid)|grepl("HS33", Plasmid), Plasmid, "")]
plasmid_tab <- plasmid_tab[,list(Plasmid, lab1, genome, Continent, Clust2)]

write.tree(plasmid_phy, file = "plasmid_cluster2/conjugative_plasmid_tree_boot.tre")

## Visualize tree with bootstraps
lab1_colors <- brewer.pal(3, "Dark2")[c(1,1,2,2,3,3)]
names(lab1_colors) <- plasmid_tab[lab1 != "", sort(Plasmid)]
plasmid_phy <- drop.tip(plasmid_phy, "MGYG000112002__contig_7") # Drop very long branch
plasmid_phy <- unroot(plasmid_phy)
plasmid_tree <- ggtree(plasmid_phy) %<+% plasmid_tab + geom_tippoint(shape = 21, aes(fill = ifelse(Clust2 == 2, "Group B", "Group A"))) + 
  geom_tiplab(size = 3, aes(label = lab1, color = lab1)) +
  scale_color_manual(values = lab1_colors, guide = "none") +
  ggnewscale::new_scale_color() +
  geom_nodelab(aes(subset = as.numeric(label) > 79, label = label), size = 3) +
  scale_fill_manual(values = brewer.pal(5, "Blues")[c(2,5)], name = "", guide = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom") + geom_treescale()
ggsave(plasmid_tree, file = "plasmid_cluster2/full_plasmid_tree_viz_boot.pdf", width = 6, height = 8)


## Other checks of plasmid seqs
consen_mat <- consensusMatrix(full_align)
conseq <- ConsensusSequence(full_align)
majority_gaps <- consen_mat[,which(apply(consen_mat, 2, function(x){
  x[names(x)=="-"] > 33
}))]
majority_gap_pos <- which(apply(consen_mat, 2, function(x){
  x[names(x)=="-"] > 33
}))
conserved_sites <- consen_mat[,which(apply(consen_mat, 2, function(x){
  x[names(x)=="-"] < 33
}))]
tot_deviations <- apply(conserved_sites, 2, function(x){
  maj_var <- x[x >= 33]
  min_vars <- x[x != 0 & x < 33]
  return(sum(min_vars))
})
1-sum(tot_deviations)/ncol(conserved_sites)/66 # ANI
has_insertion <- sapply(full_align, function(x){
  foo <- extractAt(x, IRanges(start = majority_gap_pos, width = 1))
  return(length(foo[foo != "-"]))
})
has_insertion[has_insertion > 100]
length(has_insertion[has_insertion > 100])


## Define cluster markers, cleaned up
## metagenomics/core gene snps
library(Biostrings)
library(data.table)
library(tidyverse)

# Function to ID relevant markers based on alignment
id_clust_marker_snps <- function(aln){
  consMat <- consensusMatrix(aln)
  clust2_sub <- aln[gsub(";.*", "", names(aln)) %in% metadata_qc[Completeness > 98 & Contamination < 2 & Clust2 == 2 & grepl("lenta", species), genome]]
  clust2Mat <- consensusMatrix(clust2_sub)
  clust1_sub <- aln[gsub(";.*", "", names(aln)) %in% metadata_qc[Completeness > 98 & Contamination < 2 & Clust2 == 1 & grepl("lenta", species), genome]]
  clust1Mat <- consensusMatrix(clust1_sub)
  
  # We need to identify SNPs that are only ever seen in one cluster and not the other (private)
  prev_clust1 <- t(clust1Mat[c(1:4,16),]) %>% data.table()
  prev_clust1[,Position:=1:nrow(prev_clust1)]
  prev_clust1[,Cluster:=1]
  prev_clust2 <- t(clust2Mat[c(1:4, 16),]) %>% data.table()
  prev_clust2[,Position:=1:nrow(prev_clust2)]
  prev_clust2[,Cluster:=2]
  prev_clust <- merge(prev_clust1, prev_clust2, by = "Position")
  prev_clust[,PrivateSNP:=case_when(
    A.x == 0 & A.y > 5 ~ "Clust2_A",
    A.y == 0 & A.x > 5 ~ "Clust1_A",
    C.x == 0 & C.y > 5 ~ "Clust2_C", 
    C.y == 0 & C.x > 5 ~ "Clust1_C",
    G.x == 0 & G.y > 5 ~ "Clust2_G", 
    G.y == 0 & G.x > 5 ~ "Clust1_G",
    T.x == 0 & T.y > 5 ~ "Clust2_T",
    T.y == 0 & T.x > 5 ~ "Clust1_T"
  )]
  priv_snps <- prev_clust[!is.na(PrivateSNP)]
  if(nrow(priv_snps) > 0){
    priv_snps[,c("ClusterMarker", "SNPMarker"):=tstrsplit(PrivateSNP, split = "_")]
    # Also add in prevalence info
    priv_snps[,ClustPrevalence:=sapply(1:nrow(priv_snps), function(x){
      if(priv_snps[x, ClusterMarker] == "Clust1"){
        testcol <- paste0(priv_snps[x, SNPMarker], ".x")
        totgenomes <- metadata_qc[Completeness > 98 & Contamination < 2 & grepl("lenta", species),sum(Clust2 == 1, na.rm = T)]
      } else { 
        testcol <- paste0(priv_snps[x, SNPMarker], ".y")
        totgenomes <- metadata_qc[Completeness > 98 & Contamination < 2 & grepl("lenta", species),sum(Clust2 == 2, na.rm = T)]
      }
      return(priv_snps[x, get(testcol)]/totgenomes)
    })]
  }
  
  return(priv_snps)
}


#  Import metadata
metadata_qc <- readRDS("final_qc_metadata.rds")
clust2_genomes <- metadata_qc[Clust2 == 2, genome]

## ID markers for each core gene separately
dir.create("core_gene_consensus_seqs/")
core_genes_list <- gsub(".aln.fas", "", list.files("aligned_gene_sequences/"), fixed = T)
priv_genes_all <- data.table()
genome_counts <- data.table(genome = metadata_qc[Completeness > 98 & Contamination < 2 & grepl("lenta", species), genome], 
                            Clust1Markers = 0, Clust2Markers = 0, HighPrevClust1Markers = 0, HighPrevClust2Markers = 0)

# Get markers from each core gene
for(i in 1:length(core_genes_list)){
  foo <- readDNAStringSet(paste0("aligned_gene_sequences/", core_genes_list[i], ".aln.fas"))
  foo <- foo[gsub(";.*", "", names(foo)) %in% metadata_qc[Completeness > 98 & Contamination < 2, genome]]
  consMat_f <- consensusMatrix(foo)
  cons_seq <- consensusString(consMat_f, threshold = 0.5) 
  gap_cons <- str_count(cons_seq, "-")
  if(gap_cons < 0.3*nchar(cons_seq)){ # Only proceed if reasonable number of gaps, otherwise skip this gene
    ambiguous_pos <- str_locate_all(cons_seq, "\\?")
    ambiguous_pos <- sapply(ambiguous_pos, function(x){ return(x[,"start"])})
    if(length(unlist(ambiguous_pos)) > 0){
      cons_seq <- strsplit(cons_seq, split = "")[[1]]
      for(j in 1:length(ambiguous_pos)){
        base_options <- consMat_f[,ambiguous_pos[j]]
        base_options <- names(base_options[which(base_options==max(base_options))])
        cons_seq[ambiguous_pos[j]] <- sample(base_options, 1) #Random
      }
    } else {
      cons_seq <- strsplit(cons_seq, split = "")[[1]]
    }
    #Remove gaps from both consensus and alignment
    gap_pos <- which(cons_seq == "-")
    if(length(gap_pos) > 0){
      cons_seq <- cons_seq[-gap_pos]
      foo2 <- DNAStringSet(sapply(foo, function(x){
        x2 <- as.vector(x)
        x2 <- x2[-gap_pos]
        return(paste(as.character(x2), collapse = ""))
      }))
      foo <- foo2
    }
    # Give consensus a nice name and save
    cons_seq <- DNAStringSet(paste(cons_seq, collapse = ""))
    names(cons_seq) <- paste0(core_genes_list[i], "_ElentaConsensus")
    writeXStringSet(DNAStringSet(cons_seq), file = paste0("core_gene_consensus_seqs/", 
                                                          core_genes_list[i], ".fna"))
    ## ID marker SNPs using gap-removed alignment
    priv_snps_gene <- id_clust_marker_snps(foo)
    if(nrow(priv_snps_gene) > 0){
      priv_snps_gene[,Gene:=core_genes_list[i]]
      priv_genes_all <- rbind(priv_genes_all, priv_snps_gene)
      ## Count number of markers for each genome
      for(k in 1:nrow(priv_snps_gene)){
        foo_sub <- subseq(foo, start = priv_snps_gene[k, Position], width = 1)
        foo_sub_match <- foo_sub[foo_sub == priv_snps_gene[k, SNPMarker]]
        if(priv_snps_gene[k, ClusterMarker] == "Clust1"){
          genome_counts[genome %in% gsub(";.*", "", names(foo_sub_match)), Clust1Markers:=Clust1Markers + 1]
          if(priv_snps_gene[k, ClustPrevalence] > 0.1){
            genome_counts[genome %in% gsub(";.*", "", names(foo_sub_match)), HighPrevClust1Markers:=HighPrevClust1Markers + 1]
          }
        } else{
          genome_counts[genome %in% gsub(";.*", "", names(foo_sub_match)), Clust2Markers:=Clust2Markers + 1]
          if(priv_snps_gene[k, ClustPrevalence] > 0.1){
            genome_counts[genome %in% gsub(";.*", "", names(foo_sub_match)), HighPrevClust2Markers:=HighPrevClust2Markers + 1]
          }
        }
      }
    }
  }
}
genome_counts <- merge(genome_counts, metadata_qc[,list(genome, Clust2)])
#Add ref base to priv gene table
system("cat core_gene_consensus_seqs/*.fna > core_gene_consensus_all.fna")
consSeqAll <- readDNAStringSet("core_gene_consensus_all.fna")
for(j in 1:length(consSeqAll)){
  vector_foo <- consSeqAll[[j]] %>% as.vector() %>% as.character()
  pos_keep <- priv_genes_all[Gene == gsub("_ElentaConsensus", "", names(consSeqAll)[j]), Position]
  priv_genes_all[Gene == gsub("_ElentaConsensus", "", names(consSeqAll)[j]), ConsensusRefBase:=vector_foo[pos_keep]]
}
fwrite(priv_genes_all, file = "snps_clust_markers_core_genes.tsv")
fwrite(genome_counts, file = "snps_clust_markers_counts_per_genome.tsv")


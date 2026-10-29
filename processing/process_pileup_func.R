## Final version of MGS processing function
library(data.table) 
library(tidyverse)
library(Biostrings)

process_pileup <- function(pileup_test_fil, cov_threshold = 4){
  pileup_test <- fread(pileup_test_fil)
  if(nrow(pileup_test) > 0){
    pileup_test[,Gene:=gsub("_ElentaConsensus", "", V1)]
    ## Remove insertions
    pileup_test[,SimplifiedString:=V5]
    pileup_test[grepl("+1", V5, fixed = T), SimplifiedString:=str_remove_all(V5, "\\+1[a-zA-Z]")]
    pileup_test[grepl("+2", V5, fixed = T), SimplifiedString:=str_remove_all(SimplifiedString, "\\+2[a-zA-Z][a-zA-Z]")]
    pileup_test[grepl("+3", V5, fixed = T), SimplifiedString:=str_remove_all(SimplifiedString, "\\+3[a-zA-Z][a-zA-Z][a-zA-Z]")]
    pileup_test[grepl("+4", V5, fixed = T), SimplifiedString:=str_remove_all(SimplifiedString, "\\+4[a-zA-Z]{4}")]
    pileup_test[grepl("+5", V5, fixed = T), SimplifiedString:=str_remove_all(SimplifiedString, "\\+5[a-zA-Z]{5}")]
    pileup_test[grepl("+6", V5, fixed = T), SimplifiedString:=str_remove_all(SimplifiedString, "\\+6[a-zA-Z]{6}")]
    pileup_test[grepl("+", SimplifiedString, fixed = T), SimplifiedString:=NA] ##Exclude these with long insertions
    pileup_test[,NumCharV5:=str_count(V5, pattern = "[ACGTacgt\\*]")]
    pileup_test[,NumCharSimplified:=str_count(SimplifiedString, pattern = "[ACGTacgt\\*]")]
    pileup_test[NumCharSimplified != V4 & NumCharSimplified < 8, SimplifiedString:=NA] ## Exclude those with other weird stuff/bad alignments
    pileup_test[,NumAs:=str_count(SimplifiedString, pattern = "[Aa]")]
    pileup_test[,NumCs:=str_count(SimplifiedString, pattern = "[Cc]")]
    pileup_test[,NumGs:=str_count(SimplifiedString, pattern = "[Gg]")]
    pileup_test[,NumTs:=str_count(SimplifiedString, pattern = "[Tt]")]
    pileup_test[,SimplifiedString2:=str_remove_all(toupper(SimplifiedString), pattern = "[\\+\\]\\^\\$\\!\\*\\-N0123456789]")]
    pileup_test[,NumCharSimplified:=str_count(SimplifiedString2, pattern = "[ACGTacgt]")]
    
    pileup_test[,MajorAllele:=case_when(
      NumAs >= NumCs & NumAs >= NumGs & NumAs >= NumTs ~ "A",
      NumCs >= NumAs & NumCs >= NumGs & NumCs >= NumTs ~ "C",
      NumGs >= NumCs & NumGs >= NumAs & NumGs >= NumTs ~ "G",
      NumTs >= NumCs & NumTs >= NumGs & NumTs >= NumAs ~ "T"
    )]
    pileup_test[,NumChars:=(NumAs > 0) + (NumCs > 0) + (NumGs > 0) + (NumTs > 0)]
    pileup_test[,MajorAllelePrev:=str_count(SimplifiedString2, pattern = MajorAllele)/V4]
    pileup_test[,MajorAlleleCount:=str_count(SimplifiedString2, pattern = MajorAllele)]
    pileup_test[MajorAllelePrev < 1,MinorAlleleSub:=str_remove_all(SimplifiedString2, pattern =MajorAllele)]
    pileup_test[!is.na(MinorAlleleSub) & NumChars == 2, MinorAllele:=sapply(MinorAlleleSub, function(x){
      rawToChar(unique(charToRaw(x)))
    })]
    if(nrow(pileup_test[NumChars > 2]) > 0){
      pileup_test[NumChars > 2, MinorAllele:=sapply(MinorAlleleSub, function(x){
        uniquechars <- rawToChar(unique(charToRaw(x))) %>% strsplit(split = "")
        foocount <- sapply(uniquechars[[1]], function(y){
          str_count(x, y)
        })
        return(names(foocount)[which.max(foocount)])
      })]
      pileup_test[,MinorAlleleCount:=str_count(SimplifiedString2, pattern = MinorAllele)]
      pileup_test[,MinorAllelePrev:=str_count(SimplifiedString2, pattern = MinorAllele)/V4]
    } else {
      pileup_test[,MinorAllele:=NA]
      pileup_test[,MinorAlleleCount:=NA]
      pileup_test[,MinorAllelePrev:=NA]
    }
    pileup_test[is.na(MinorAllele) & !is.na(MajorAllele), MinorAlleleCount:=0]
    pileup_test[is.na(MinorAllele)& !is.na(MajorAllele), MinorAllelePrev:=0]
    setnames(pileup_test, "V2", "Position")
    pileup_test_priv <- merge(pileup_test, priv_genes_all, by = c("Gene", "Position"), all.x = F, all.y = T) #Missing ones will be NA but the filter on the next step still works OK
    ## Save only high-prevalence markers (10% min)
    pileup_test_conf <- pileup_test_priv[NumCharSimplified > cov_threshold & ClustPrevalence > 0.1]
    pileup_test_conf[,Sample:=gsub("_pileup_out.txt", "", basename(pileup_test_fil))]
    return(pileup_test_conf)
  }
}
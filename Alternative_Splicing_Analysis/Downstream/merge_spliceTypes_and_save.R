merge_spliceTypes_and_save <- function(splice_events, new_suffix,indir){ 
  
  ##""" Take modified tables with Replicates-columns & Exon-coordinates 
  ##""" and harmonize columns for all splicing-events
  
  listSpliceTabs <- list()
  for (s in splice_events) {
    splice_tab <- read.delim(paste0(indir, s, new_suffix),header=T,as.is=T)
    colsOI <- c("GeneID",	"geneSymbol","chr","strand",
                "IncFormLen", "SkipFormLen","ExonLength",
                grep("IncLevel1_*", colnames(splice_tab),value=T),
                grep("IncLevel2_*", colnames(splice_tab),value=T),
                "splice_event", "IncLevelDifference",
                "alternative_exon_coordinates","PValue","FDR" )
    
    listSpliceTabs[[s]] <-  splice_tab[,colsOI]
  }
  
  ##---SAVE ----##
  merged_splice_events <- do.call(rbind, listSpliceTabs)
 
  return(merged_splice_events)
  
}
#!/usr/bin/env R 

###> Author: Maria Kondili
###> Subject :Merge using all records of each ase-table,not only common exon-coord.
###> 


suppressPackageStartupMessages(library(tidyverse))

full_ase <- dplyr::full_join(ase_mrg_NT, 
                             ase_mrg_sg2,
                             by="alternative_exon_coordinates",
                             suffix=c(".x",".y"))

full_ase <- full_ase %>%  dplyr::rename("dPSI_CTvsDM1NT"    = "dPSI.x",
                                        "dPSI_DM1_sgDMPK"   = "dPSI.y",
                                        "PValue_CTvsDM1NT"  = "PValue.x",
                                        "PValue_DM1_sgDMPK" = "PValue.y",
                                        "FDR_CTvsDM1NT"     = "FDR.x",
                                        "FDR_DM1_sgDMPK"    = "FDR.y" )
##> all_ase : 201
##> full_ase: 106773 ( sans filtrage FDR,dPSI )


## see which events are added from y that are not common to x-obj :
full_ase$GeneID.y %>% is.na %>% which %>% length # 13012 

full_ase <- dplyr::full_join(full_ase, ase_mrg_ct_vs_sg2, 
                             by="alternative_exon_coordinates", 
                             suffix = c(".z",".w"))

full_ase <- full_ase %>%  dplyr::rename("dPSI_CTvsDMsgDMPK" = "dPSI",
                                        "PValue_CTvsDMsgDMPK" = "PValue",
                                        "FDR_CTvsDMsgDMPK"="FDR" )

###>Add Correction, calculated from other script in this table
###> all_ase$Correction_Pct

all_ase <- read_delim("Results/Merged_common_ASevents_with_CorrectionRatio_CTnt_DM1nt_DM1sgDMPK.tsv",
            delim="\t",col_names=T)

ase_with_corx <- match(full_ase$alternative_exon_coordinates ,all_ase$alternative_exon_coordinates)
names(ase_with_corx) <- seq(1,nrow(full_ase))
ase_with_corx <- ase_with_corx[-which(is.na(ase_with_corx) )]

#initialise Correction:
#full_ase <- mutate(full_ase, "Correction_Pct"= NA)
#full_ase <- mutate(full_ase, "CorrectionRatio" = NA)

# add values to column "Correction"
full_ase[as.integer(names(ase_with_corx)),"Correction_Ratio"] <- all_ase[ase_with_corx,"Correction_Ratio"]

write_delim(x=full_ase,
            file=paste0(final_outdir_all,"Merged_all_AS_events_withCorrection_full_join.tsv"),
            delim="\t",col_names = TRUE)









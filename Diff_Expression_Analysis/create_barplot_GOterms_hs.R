create_barplot_GOterms.hs <- function(genenames,gene_format = "SYMBOL",
                                      my_ontol="ALL",plots_dir,main="") {

  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  suppressPackageStartupMessages(library(AnnotationDbi))

  EntrezID <- mapIds(x=org.Hs.eg.db, keys=genenames,
                     column="ENTREZID",
                     keytype=gene_format,
                     multiVals = "first")
  ##to see accepted keytypes :  keytypes(org.Mm.eg.db)
  ## column = what I want to obtain
  ## keytype = what I give in keys

  suppressPackageStartupMessages(library(DOSE))
  suppressPackageStartupMessages(library(enrichplot))
  go.enrichm.bp <- enrichGO(EntrezID, ont=my_ontol,
                            OrgDb=org.Hs.eg.db,
                            keyType = "ENTREZID",
                            pvalueCutoff=0.1,
                            qvalueCutoff=0.1,
                            pAdjustMethod="BH",
                            readable=TRUE)

  # go.enrichm.mf <- enrichGO(sig.data.d$Entrez_ID, ont="MF",OrgDb=org.Mm.eg.db,keyType = "ENTREZID",
  #                           pvalueCutoff=0.05, qvalueCutoff=0.1, pAdjustMethod = "BH", readable=TRUE)
  #

  if (my_ontol=="BP") {msg = "Biol_Processes" }
  if (my_ontol=="MF") {msg = "Molec_Functions"}
  if (my_ontol=="CC") {msg = "Cellular_Components" }
  if (my_ontol== "ALL") {msg="All_GO_terms"}

   png(paste0(plots_dir,"Barplot_",msg,"_",main,".png"))
    p <- barplot(go.enrichm.bp, showCategory = 20,
            x="GeneRatio", color="p.adjust", font.size=8,
            title=paste(msg," of ", main))
    #p
    print(p)
   dev.off()
}

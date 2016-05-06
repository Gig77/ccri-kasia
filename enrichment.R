library(GO.db)
library(clusterProfiler)

enrichAll <- function(gene.entrez=NULL, gene.hgnc=NULL, universe.entrez=NULL, universe.hgnc=NULL, pvalueCutoff = 0.05, qvalueCutoff = 0.2) {
  
  result <- data.frame(Category=character(0), ID=character(0), Description=character(0), GeneRatio=character(0), BgRatio=character(0), pvalue=numeric(0), p.adjust=numeric(0), qvalue=numeric(0), Count=integer(0))
  
  universe.mf <- unique(org.Hs.egGO2ALLEGS[["GO:0003674"]])
  ego.mf <- enrichGO(
    gene = gene.entrez,
    universe = if(is.null(universe.entrez)) universe.mf else unique(intersect(universe.entrez, universe.mf)),
    organism = "human",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    readable = TRUE
  )
  if (nrow(ego.mf@result) > 0) result <- rbind(result, cbind(data.frame(Category="GO molecular function"), ego.mf@result))
  
  universe.bp <- unique(org.Hs.egGO2ALLEGS[["GO:0008150"]])
  ego.bp <- enrichGO(
    gene = gene.entrez,
    universe = if(is.null(universe.entrez)) universe.bp else unique(intersect(universe.entrez, universe.bp)),
    organism = "human",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    readable = TRUE
  )
  if (nrow(ego.bp@result) > 0) result <- rbind(result, cbind(data.frame(Category="GO biological process"), ego.bp@result))
  
  universe.cc <- unique(org.Hs.egGO2ALLEGS[["GO:0005575"]])
  ego.cc <- enrichGO(
    gene = gene.entrez,
    universe = if(is.null(universe.entrez)) universe.cc else unique(intersect(universe.entrez, universe.cc)),
    organism = "human",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    readable = TRUE
  )
  if (nrow(ego.cc@result) > 0) result <- rbind(result, cbind(data.frame(Category="GO cellular component"), ego.cc@result))
  
  kegg <- enrichKEGG(
    gene = gene.entrez,
    universe = universe.entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    readable = TRUE
  )
  if (nrow(kegg@result) > 0) result <- rbind(result, cbind(data.frame(Category="KEGG"), kegg@result))
  
  msigdb <- enricher(  
    gene          = gene.hgnc, 
    pAdjustMethod = "BH", 
    universe      = universe.hgnc, 
    minGSSize     = 5, 
    pvalueCutoff  = pvalueCutoff,
    qvalueCutoff  = qvalueCutoff,
    TERM2GENE     = term2gene, 
    TERM2NAME     = gmt[,c(1,2)]
  )
  if (nrow(msigdb@result) > 0) result <- rbind(result, cbind(data.frame(Category="MSigDB5.0"), msigdb@result))
  
  result[order(result$pvalue),]
}

# MSigDB 5.0 gene sets
gmt <- read.delim("/mnt/projects/generic/data/msigdb5.0/msigdb.v5.0.symbols.gmt", stringsAsFactors = F, header = FALSE)
term2gene <- do.call("rbind", apply(gmt[gmt$V3 != "",], 1, function(x) {
  data.frame(term=x[1], gene=unique(x[3:length(x)][x[3:length(x)]!=""]), row.names = NULL)
}))

# get background entrez gene ids (=genes expressed in both fuka and kasia dataset)

library(biomaRt)
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice" , dataset="hsapiens_gene_ensembl") # GRCh37, v75
fuka.d20plus.REH <- read.delim("/mnt/projects/chrisi/results/fuka/telamlKD.REHandAT2.esetnsF.onlyG_late_REH.annot.tsv")
kasia.H24vsC24 <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_H24vsC24/table.csv")
kasia.H48vsC48 <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_H48vsC48/table.csv")
bg.hgnc <- unique(intersect(fuka.d20plus.REH$syms, union(kasia.H24vsC24$Gene[!is.na(kasia.H24vsC24$q)], kasia.H48vsC48$Gene[!is.na(kasia.H48vsC48$q)])))
bg.entrez <- as.character(unique(getBM(attributes=c("entrezgene"), filters = "hgnc_symbol", values=bg.hgnc, mart=mart)$entrezgene))

# REH ChIP-seq peaks

reh <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)

# filter for fuka up-regulated and promoter peaks

reh.FUKAup <- with(reh, reh[!is.na(fuka.kdER.REH.d20plus.padj) & fuka.kdER.REH.d20plus.padj <= 0.2 & !is.na(fuka.kdER.REH.d20plus.logfc) & fuka.kdER.REH.d20plus.logfc > 0,])
reh.FUKAup.prom <- with(reh.FUKAup, reh.FUKAup[`Distance to TSS` >= -5000 & `Distance to TSS` <= 2000,])

# H24

reh.FUKAup.prom.H24up <- with(reh.FUKAup.prom, reh.FUKAup.prom[!is.na(kasia.HDACi.H24vsC24.padj) & kasia.HDACi.H24vsC24.padj <= 0.2 & !is.na(kasia.HDACi.H24vsC24.logfc) & kasia.HDACi.H24vsC24.logfc > 0,])

reh.FUKup.prom.H24up.enr <- enrichAll(
  gene.entrez=as.character(unique(reh.FUKAup.prom.H24up$`Entrez ID`)),
  gene.hgnc=unique(reh.FUKAup.prom.H24up$`Gene Name`),
  universe.entrez=bg.entrez,
  universe.hgnc=bg.hgnc,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

write.table(reh.FUKAup.prom.H24up, "/mnt/projects/kasia/results/genes.H24-up.fuka-up.promoter-peak.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(reh.FUKup.prom.H24up.enr, "/mnt/projects/kasia/results/enrichment.H24-up.fuka-up.promoter-peak.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

# H48

reh.FUKAup.prom.H48up <- with(reh.FUKAup.prom, reh.FUKAup.prom[!is.na(kasia.HDACi.H48vsC48.padj) & kasia.HDACi.H48vsC48.padj <= 0.2 & !is.na(kasia.HDACi.H48vsC48.logfc) & kasia.HDACi.H48vsC48.logfc > 0,])

reh.FUKup.prom.H48up.enr <- enrichAll(
  gene.entrez=as.character(unique(reh.FUKAup.prom.H48up$`Entrez ID`)),
  gene.hgnc=unique(reh.FUKAup.prom.H48up$`Gene Name`),
  universe.entrez=bg.entrez,
  universe.hgnc=bg.hgnc,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

write.table(reh.FUKAup.prom.H48up, "/mnt/projects/kasia/results/genes.H48-up.fuka-up.promoter-peak.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(reh.FUKup.prom.H48up.enr, "/mnt/projects/kasia/results/enrichment.H48-up.fuka-up.promoter-peak.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

kasia.H24vsC24 <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_H24vsC24/table.csv", stringsAsFactors = F)
kasia.H24vsC24 <- kasia.H24vsC24[order(kasia.H24vsC24$q),]
kasia.H24vsC24 <- kasia.H24vsC24[!duplicated(kasia.H24vsC24$Gene),]
kasia.H24vsC24 <- kasia.H24vsC24[,c("Gene", "q", "fc")]
names(kasia.H24vsC24) <- c("Gene", "kasia.HDACi.H24vsC24.padj", "kasia.HDACi.H24vsC24.logfc")
combined <- kasia.H24vsC24

kasia.H48vsC48 <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_H48vsC48/table.csv", stringsAsFactors = F)
kasia.H48vsC48 <- kasia.H48vsC48[order(kasia.H48vsC48$q),]
kasia.H48vsC48 <- kasia.H48vsC48[!duplicated(kasia.H48vsC48$Gene),]
kasia.H48vsC48 <- kasia.H48vsC48[,c("Gene", "q", "fc")]
names(kasia.H48vsC48) <- c("Gene", "kasia.HDACi.H48vsC48.padj", "kasia.HDACi.H48vsC48.logfc")
combined <- merge(combined, kasia.H48vsC48, all.x = T)

fiona.peaks.REH <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_peaks.annotated.tsv", check.names = F, stringsAsFactors = F)
fiona.peaks.REH <- fiona.peaks.REH[fiona.peaks.REH$`Distance to TSS` >= -5000 & fiona.peaks.REH$`Distance to TSS` <= 2000, c("Gene Name", "Peak Score")]
fiona.peaks.REH <- fiona.peaks.REH[order(fiona.peaks.REH$`Peak Score`, decreasing = T),]
fiona.peaks.REH <- fiona.peaks.REH[!duplicated(fiona.peaks.REH$`Gene Name`),]
names(fiona.peaks.REH) <- c("Gene", "fiona.REH.promoterpeak.score")
combined <- merge(combined, fiona.peaks.REH, all.x = T)

fuka.d20plus.REH <- read.delim("/mnt/projects/chrisi/results/fuka/telamlKD.REHandAT2.esetnsF.onlyG_late_REH.annot.tsv")
fuka.d20plus.REH <- fuka.d20plus.REH[,c("syms", "Padj.onlyG_late_REH", "logFC.onlyG_late_REH")]
names(fuka.d20plus.REH) <- c("Gene", "fuka.kdER.d20plus.REH.padj", "fuka.kdER.d20plus.REH.logfc")
combined <- merge(combined, fuka.d20plus.REH, all.x = T)

# plot H24

pdf("/mnt/projects/kasia/results/kasia-H24-vs-Fuka-REH-D20+.pdf")

sigLevel <- 0.2
combined.sig <- combined[!is.na(combined$kasia.HDACi.H24vsC24.padj) & combined$kasia.HDACi.H24vsC24.padj <= sigLevel & !is.na(combined$fuka.kdER.d20plus.REH.logfc),]
fit <- lm(fuka.kdER.d20plus.REH.logfc~kasia.HDACi.H24vsC24.logfc, combined.sig)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(fuka.kdER.d20plus.REH.logfc~kasia.HDACi.H24vsC24.logfc, 
     data=combined.sig, 
     xlab="log2FC HDACi H24", 
     ylab="log2FC Fuka REH D20+", 
     ylim=c(-3, 3), 
     xlim=c(-3, 3), 
     main=sprintf("HDACi H24 q<=%.2g (R=%.2f, p=%.2g)", sigLevel, R, p), 
     cex=0.4,
     col=ifelse(is.na(combined.sig$fiona.REH.promoterpeak.score), "black", "red"))
abline(fit, col="red")
abline(v=0, lty=3)
abline(h=0, lty=3)
combined.sig$label <- with(combined.sig, Gene %in% unique(c(Gene[order(kasia.HDACi.H24vsC24.logfc, decreasing=F)][1:10],
                                                            Gene[order(fuka.kdER.d20plus.REH.logfc, decreasing=F)][1:10],
                                                            Gene[order(kasia.HDACi.H24vsC24.logfc, decreasing=T)][1:10],
                                                            Gene[order(fuka.kdER.d20plus.REH.logfc, decreasing=T)][1:10])))
with(combined.sig[combined.sig$label,], 
     text(kasia.HDACi.H24vsC24.logfc,
          fuka.kdER.d20plus.REH.logfc + 0.07, 
          Gene,
          cex=0.4))

legend("bottomright", c("w/o E/R promoter peak in REH", "w/ E/R promoter peak in REH"), fill=c("black", "red"), cex=0.8)

dev.off()

# plot H48

pdf("/mnt/projects/kasia/results/kasia-H48-vs-Fuka-REH-D20+.pdf")

sigLevel <- 0.2
combined.sig <- combined[!is.na(combined$kasia.HDACi.H48vsC48.padj) & combined$kasia.HDACi.H48vsC48.padj <= sigLevel & !is.na(combined$fuka.kdER.d20plus.REH.logfc),]
fit <- lm(fuka.kdER.d20plus.REH.logfc~kasia.HDACi.H48vsC48.logfc, combined.sig)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(fuka.kdER.d20plus.REH.logfc~kasia.HDACi.H48vsC48.logfc, 
     data=combined.sig, 
     xlab="log2FC HDACi H48", 
     ylab="log2FC Fuka REH D20+", 
     ylim=c(-3, 3), 
     xlim=c(-3, 3), 
     main=sprintf("HDACi H48 q<=%.2g (R=%.2f, p=%.2g)", sigLevel, R, p), 
     cex=0.4,
     col=ifelse(is.na(combined.sig$fiona.REH.promoterpeak.score), "black", "red"))
abline(fit, col="red")
abline(v=0, lty=3)
abline(h=0, lty=3)
combined.sig$label <- with(combined.sig, Gene %in% unique(c(Gene[order(kasia.HDACi.H48vsC48.logfc, decreasing=F)][1:10],
                                                            Gene[order(fuka.kdER.d20plus.REH.logfc, decreasing=F)][1:10],
                                                            Gene[order(kasia.HDACi.H48vsC48.logfc, decreasing=T)][1:10],
                                                            Gene[order(fuka.kdER.d20plus.REH.logfc, decreasing=T)][1:10])))
with(combined.sig[combined.sig$label,], 
     text(kasia.HDACi.H48vsC48.logfc,
          fuka.kdER.d20plus.REH.logfc + 0.07, 
          Gene,
          cex=0.4))

legend("bottomright", c("w/o E/R promoter peak in REH", "w/ E/R promoter peak in REH"), fill=c("black", "red"), cex=0.8)

dev.off()
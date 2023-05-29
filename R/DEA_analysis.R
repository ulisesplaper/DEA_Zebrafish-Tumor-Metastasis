# Cargar librerias necesarias
library(limma)
library(edgeR)
library("pheatmap")

#Importar el archivo de fc
fc <- read.table("data/counts.txt", header = T)

# Cargar informacion al objeto de edgeR
samples <- factor(c("Parental", "Parental", "Metastatic", "Metastatic"))
DGEList_ZTM = DGEList(counts=fc[,3:6],
                      group=samples,
                      genes=fc[,1:2])


# Eliminar genes con una baja expresion
keep = rowSums(cpm(DGEList_ZTM)>1) >= 3
DGEList_ZTM = DGEList_ZTM[keep,]

# Generar matriz de diseño
design <- model.matrix(~0+samples)
colnames(design) <- levels(samples)

# Factores de normalizacion y variacion de la muestra
DGEList_ZTM = calcNormFactors(DGEList_ZTM)
DGEList_ZTM = estimateDisp(DGEList_ZTM, design=design)

colnames(DGEList_ZTM$counts) <- c("P1","P2", "M1","M2")

# Generar un pca
color_all = c("blue","blue","red","red")
plotMDS(DGEList_ZTM, cex=0.8,col=color_all, )
dev.copy(pdf,"MDS-ZTM_RNA.pdf")
dev.off()

# Calcular genes diferencialmente expresados
fit <- glmFit(DGEList_ZTM, design)
contrasts = makeContrasts(ParentalvsMetastatic=Metastatic-Parental, levels=design)
lrt_hip <- glmLRT(fit, contrast=contrasts[,"ParentalvsMetastatic"])
DEG_lrt <- as.data.frame(topTags(lrt_hip, n=length(DGEList_ZTM$counts[,1])))

# Extraer genes upregulated
up = (DEG_lrt$logFC > 2) & (DEG_lrt$FDR < 1e-5)
write(DEG_lrt[up,1], "upregulategenes_ZTM.txt")
# Extraer los IDs de los genes upregulated
write.table(DEG_lrt[up,], "edgeR-DEG-MetastaticvsParental-Up--1e-5.txt", sep="\t", quote=FALSE,row.names=FALSE)

# Extraer genes downregulated
down = (DEG_lrt$logFC < -2) & (DEG_lrt$FDR < 1e-5)
write(DEG_lrt[down,1], "downregulategenes_ZTM.txt")
# Extrae los IDs de los genes upregulated
write.table(DEG_lrt[down,], "edgeR-DEG-MetastaticvsParental-Down-1e-5.txt", sep="\t", quote=FALSE,
            row.names=FALSE)

# Volcano plot
plot(DEG_lrt$logFC[!(up | down)], -
       log10(DEG_lrt$FDR[!(up | down)]), pch=19,
     col="gray", cex=0.4, xlab="log2 Expression fold
change", ylab="-log10 FDR", main="Volcano plot
Parental vs Metastatic", xlim=c(-10,10),ylim=c(0,30))
points(DEG_lrt$logFC[up], -
         log10(DEG_lrt$FDR[up]), pch=19, col="red",
       cex=0.4)
points(DEG_lrt$logFC[down], -
         log10(DEG_lrt$FDR[down]), pch=19, col="blue",
       cex=0.4)
abline(h=5, col="black", lty=3)
abline(v=c(-2,2), col="black", lty=3)
dev.copy(pdf, "Volcano_plot-RNA_seqWTvsKO.pdf")
dev.off()


log2_fpkm_ZTM = rpkm(DGEList_ZTM,DGEList_ZTM$genes$Length, log=T)
log2_fpkm_ZTM_average = rpkmByGroup(DGEList_ZTM, log=T)

# Funcion para calcular el log PKM a partir del log fpkm
fpkm2tpm_log2 <- function(fpkm) { fpkm -log2(sum(2^fpkm)) + log2(1e6) }
log2_tpm_ZTM = apply(log2_fpkm_ZTM, 2,fpkm2tpm_log2)
log2_tpm_ZTM_average = apply(log2_fpkm_ZTM_average, 2, fpkm2tpm_log2)
# Guardar los archivos
write.table(log2_tpm_ZTM, "mm10-RNA-ZTM_log2.txt", sep="\t", quote=FALSE)
write.table(log2_tpm_ZTM_average, "mm10-RNA-ZTM-average.txt", sep="\t",quote=FALSE)

# Generar un heatmap
# Concatenar el vector de genes up y down regulated
heatmapVect <- c(rownames(DEG_lrt[up,]), rownames(DEG_lrt[down,]))

# Dataframe con los tpm
exprcol <- log2_tpm_ZTM[heatmapVect,]
upreg <- data.frame(exprcol, row.names = c(DEG_lrt[up,1], DEG_lrt[down,1]))

# Graficar el heatmap
pheatmap(
  upreg,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
    show_colnames = TRUE,legend = T
)
dev.copy(pdf, "heatmap_PvsM.pdf")
dev.off()



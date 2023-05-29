# Generar la tabla de fc (se hace en R)
library(Rsubread)
bfiles <- c("SRR12082751_P1.bam","SRR12082752_P.bam","SRR12082753_c.bam","SRR12082754_c.bam" )
fc <- featureCounts(files=bfiles,
                    annot.ext="hg38.gtf",
                    isGTFAnnotationFile=T, useMetaFeatures=T,
                    minMQS=10, largestOverlap=T, isPairedEnd=T,
                    requireBothEndsMapped=T, nthreads=5)

# Exportar la tabla para trabajar en la computadora
write.table(
  x=data.frame(fc$annotation[,c("GeneID","Length")],
               fc$counts,
               stringsAsFactors=FALSE),
  file="counts.txt",
  quote=FALSE,
  sep="\t",
  row.names=FALSE)

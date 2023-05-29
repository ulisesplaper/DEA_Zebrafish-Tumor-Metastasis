# DEA_Zebrafish-Tumor-Metastasis

## Fundamento biológico

El modelo pez cebra-tumor-metástasis (ZTM) es un método robusto para el aislamiento de células metastásicas y está especialmente indicado para facilitar el estudio de la metástasis del cáncer y las terapias de prevención de la metástasis.

## Overall design:

Poblaciones celulares parentales son extraidas del tumor primario (2 replicas). Poblaciones metastáticas son extraídas de regiones distantes al tumor primario

![Diseño experimental](https://github.com/ulisesplaper/DEA_Zebrafish-Tumor-Metastasis/blob/master/data/generaldesign.png?raw=true)

## Script summary

### 1. bash/preprocesamiento.jdl

Descarga de los archivos FASTQ, analisis de calidad, trimming, alieamiento al genoma de referencia y generación de archivos bam.

### 2. R/counts_generation.R

Genera archivos de cuentas a partir de los archivos bam.

### 3. R/DEA_analysis.R

Realiza analisis de expresion diferencial y genera mapas de interes.

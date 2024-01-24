# Carga de las librerías pertinentes.
library(clusterProfiler)
library(ggplot2)

library(MAGeCKFlute)

# Directorio de trabajo (Cambiar si fuera necesario).
setwd("Directorio/Analysis_Results")

### PatuT_Control

# Carga del archivo "all.countsummary.txt."
countsummary <- read.table(
  "PatuT_rra_control/results/count/all.countsummary.txt", 
                        header = TRUE, sep = "\t")
head(countsummary)

# Creación del directorio para los resultados de visualización de qc.
dir.create("PatuT_rra_control/results/qc/MAGeCKFlute")

# Elaboración del gráfico con los índices de Gini por muestra.
pdf("PatuT_rra_control/results/qc/MAGeCKFlute/Gini_Index.pdf")
BarView(countsummary, x = "Label", y = "GiniIndex",
        ylab = "Gini index", main = "Evenness of sgRNA reads") +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 
dev.off()

# Elaboración del gráfico con la proporción de sgRNAs sin lecturas.
pdf("PatuT_rra_control/results/qc/MAGeCKFlute/Missed_sgRNAs.pdf")
countsummary$Missed = log10(countsummary$Zerocounts)
BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "Log10 missed gRNAs", main = "Missed sgRNAs") +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 
dev.off()

# Elaboración del gráfico con el total de lecturas mapeadas.
pdf("PatuT_rra_control/results/qc/MAGeCKFlute/Mapping_Ratio.pdf")
MapRatesView(countsummary) +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()

### PatuS_Control

countsummary <- read.table(
  "PatuS_rra_control/results/count/all.countsummary.txt", 
  header = TRUE, sep = "\t")
head(countsummary)

dir.create("PatuS_rra_control/results/qc/MAGeCKFlute")

pdf("PatuS_rra_control/results/qc/MAGeCKFlute/Gini_Index.pdf")
BarView(countsummary, x = "Label", y = "GiniIndex",
        ylab = "Gini index", main = "Evenness of sgRNA reads") +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 
dev.off()

pdf("PatuS_rra_control/results/qc/MAGeCKFlute/Missed_sgRNAs.pdf")
countsummary$Missed = log10(countsummary$Zerocounts)
BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "Log10 missed gRNAs", main = "Missed sgRNAs") +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 
dev.off()

pdf("PatuS_rra_control/results/qc/MAGeCKFlute/Mapping_Ratio.pdf")
MapRatesView(countsummary) +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()

### PatuT_Day0

countsummary <- read.table(
  "PatuT_rra_day0/results/count/all.countsummary.txt", 
  header = TRUE, sep = "\t")
head(countsummary)

dir.create("PatuT_rra_day0/results/qc/MAGeCKFlute")

pdf("PatuT_rra_day0/results/qc/MAGeCKFlute/Gini_Index.pdf")
BarView(countsummary, x = "Label", y = "GiniIndex",
        ylab = "Gini index", main = "Evenness of sgRNA reads") +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 
dev.off()

pdf("PatuT_rra_day0/results/qc/MAGeCKFlute/Missed_sgRNAs.pdf")
countsummary$Missed = log10(countsummary$Zerocounts)
BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "Log10 missed gRNAs", main = "Missed sgRNAs") +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 
dev.off()

pdf("PatuT_rra_day0/results/qc/MAGeCKFlute/Mapping_Ratio.pdf")
MapRatesView(countsummary) +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()

### PatuS_Day0

countsummary <- read.table(
  "PatuS_rra_day0/results/count/all.countsummary.txt", 
  header = TRUE, sep = "\t")
head(countsummary)

dir.create("PatuS_rra_day0/results/qc/MAGeCKFlute")

pdf("PatuS_rra_day0/results/qc/MAGeCKFlute/Gini_Index.pdf")
BarView(countsummary, x = "Label", y = "GiniIndex",
        ylab = "Gini index", main = "Evenness of sgRNA reads") +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 
dev.off()

pdf("PatuS_rra_day0/results/qc/MAGeCKFlute/Missed_sgRNAs.pdf")
countsummary$Missed = log10(countsummary$Zerocounts)
BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "Log10 missed gRNAs", main = "Missed sgRNAs") +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 
dev.off()

pdf("PatuS_rra_day0/results/qc/MAGeCKFlute/Mapping_Ratio.pdf")
MapRatesView(countsummary) +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()
# 1. Importación de featureCounts: Cargar los conteos de genes obtenidos mediante alineamiento al genoma.
# Preparar una matriz con las tres muestras para comparar con Salmon.

fc <- read.delim("results/alignment_based/featurecounts/gene_counts_tail1.txt", comment.char = "#", check.names = FALSE)

fc_counts <- fc[, 7:9]
rownames(fc_counts) <- fc$Geneid
colnames(fc_counts) <- c("SRR19020908", "SRR19020909", "SRR19020910")
fc_counts <- as.matrix(fc_counts)

# 2. Detección de genes: Contar los genes con al menos una lectura asignada en cada muestra.
# Evaluar cuántos genes se detectan en al menos una muestra y en las tres muestras.

genes_detectados <- colSums(fc_counts > 0)
genes_detectados
sum(rowSums(fc_counts > 0) > 0)
sum(rowSums(fc_counts > 0) == 3)

saveRDS(fc_counts, "results/gene_expression/featurecounts_gene_counts.rds")

# 3. Importación de Salmon y tximport: Cargar los conteos estimados por Salmon y resumidos a nivel de gen mediante tximport.
# Evaluar cuántos genes se detectan en una, todas o ninguna de las tres muestras.

txi <- readRDS("results/gene_expression/tximport_gene_counts.rds")

salmon_counts <- txi$counts
dim(salmon_counts)

sum(rowSums(salmon_counts > 0) > 0)
sum(rowSums(salmon_counts > 0) == 3)
sum(rowSums(salmon_counts > 0) == 0)

# 4.  Comparación de genes entre featureCounts y Salmon: Identificar los genes presentes en ambos métodos y los exclusivos de cada uno.
# Comparar cuántos genes se detectan en ambos métodos, solo en uno o en ninguno, utilizando los umbrales > 0 y >= 1.

genes_fc <- rownames(fc_counts)
genes_salmon <- rownames(salmon_counts)

genes_comunes <- intersect(genes_fc, genes_salmon)

genes_solo_fc <- setdiff(genes_fc, genes_salmon)
genes_solo_salmon <- setdiff(genes_salmon, genes_fc)

length(genes_comunes)
length(genes_solo_fc)
length(genes_solo_salmon)

fc_comunes <- fc_counts[genes_comunes, , drop = FALSE]
salmon_comunes <- salmon_counts[genes_comunes, , drop = FALSE]

stopifnot(identical(rownames(fc_comunes), rownames(salmon_comunes)))
stopifnot(identical(colnames(fc_comunes), colnames(salmon_comunes)))

detectados_fc <- rowSums(fc_comunes > 0) > 0
detectados_salmon <- rowSums(salmon_comunes > 0) > 0

sum(detectados_fc & detectados_salmon)
sum(detectados_fc & !detectados_salmon)
sum(!detectados_fc & detectados_salmon)
sum(!detectados_fc & !detectados_salmon)

detectados_fc_1 <- rowSums(fc_comunes >= 1) > 0
detectados_salmon_1 <- rowSums(salmon_comunes >= 1) > 0

sum(detectados_fc_1 & detectados_salmon_1)
sum(detectados_fc_1 & !detectados_salmon_1)
sum(!detectados_fc_1 & detectados_salmon_1)
sum(!detectados_fc_1 & !detectados_salmon_1)


# 5. Correlación de Spearman: Evaluar la concordancia de la cuantificación entre featureCounts y Salmon para cada muestra.
# Calcular la correlación utilizando únicamente genes con conteos positivos en ambos métodos.

resultados_correlacion <- data.frame(muestra = colnames(fc_comunes), genes_comunes_detectados = NA_integer_, spearman = NA_real_)

for (i in seq_along(colnames(fc_comunes))) {
  
  fc_muestra <- fc_comunes[, i]
  salmon_muestra <- salmon_comunes[, i]
  
  positivos <- fc_muestra > 0 & salmon_muestra > 0
  
  resultados_correlacion$genes_comunes_detectados[i] <- sum(positivos)
  
  resultados_correlacion$spearman[i] <- if (sum(positivos) >= 2) cor(fc_muestra[positivos], salmon_muestra[positivos], method = "spearman") else NA_real_
}

resultados_correlacion


# 6. Visualización de la correlación: Crear tres gráficos de dispersión para comparar la cuantificación de los genes en cada muestra.
# Mostrar los genes positivos en ambos métodos, la correlación de Spearman y la línea de igualdad.

muestras <- colnames(fc_comunes)
colores <- c("#2563EB", "#E87924", "#159A75")

limite_max <- max(log1p(fc_comunes), log1p(salmon_comunes))

dir.create("results/gene_expression/figures", recursive = TRUE, showWarnings = FALSE)

png("results/gene_expression/figures/featurecounts_vs_salmon_3panels.png", width = 2400, height = 850, res = 200)

par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

for (i in seq_along(muestras)) {
  
  fc_valores <- fc_comunes[, i]
  salmon_valores <- salmon_comunes[, i]
  
  positivos <- fc_valores > 0 & salmon_valores > 0
  
  plot(log1p(fc_valores[positivos]),
       log1p(salmon_valores[positivos]),
       xlim = c(0, limite_max),
       ylim = c(0, limite_max),
       pch = 16, cex = 0.3,
       col = adjustcolor(colores[i], alpha.f = 0.2),
       xlab = "featureCounts: log(1 + counts)",
       ylab = "Salmon: log(1 + estimated counts)",
       main = paste0(muestras[i], "\nSpearman = ",
                     round(resultados_correlacion$spearman[i], 3)))
  
  abline(a = 0, b = 1, lty = 2, lwd = 1.5)
}

dev.off()


# 7. Identificación de genes discrepantes: Calcular las diferencias de cuantificación entre Salmon y featureCounts para cada gen.
# Seleccionar los cinco genes con mayor diferencia absoluta media de log2 entre las tres muestras.

diferencias_log2 <- log2((salmon_comunes + 1) / (fc_comunes + 1))
diferencia_media <- rowMeans(abs(diferencias_log2))
expresion_media <- rowMeans((salmon_comunes + fc_comunes) / 2)

comparacion_genes <- data.frame(Geneid = rownames(fc_comunes), diferencia_media_absoluta = diferencia_media, expresion_media = expresion_media)
comparacion_genes <- comparacion_genes[order(-comparacion_genes$diferencia_media_absoluta), ]

head(comparacion_genes, 10)

top5_genes <- head(comparacion_genes$Geneid, 5)

fc_comunes[top5_genes, ]
salmon_comunes[top5_genes, ]
round(diferencias_log2[top5_genes, ], 2)


# 8. Mapa de calor de genes discrepantes: Visualizar las diferencias log2 de los cinco genes con mayores discrepancias.
# Los valores positivos indican mayor cuantificación por Salmon y los negativos mayor cuantificación por featureCounts.

matriz_top5 <- diferencias_log2[top5_genes, , drop = FALSE]

colores_heatmap <- colorRampPalette(c("#3979B6", "white", "#BA5545"))(101)

nombres_genes <- c(
  "ENSG00000254772" = "EEF1G",
  "ENSG00000225840" = "processed pseudogene",
  "ENSG00000167996" = "FTH1",
  "ENSG00000232573" = "RPL3P4",
  "ENSG00000184584" = "STING1"
)

nombres_top5 <- nombres_genes[top5_genes]
nombres_top5[is.na(nombres_top5)] <- "Gene ID"

etiquetas_genes <- paste0(rownames(matriz_top5), "\n", "(", nombres_top5, ")")

png("results/gene_expression/figures/top5_gene_discrepancies_heatmap.png", width = 1800, height = 1100, res = 200)

par(mar = c(13, 12, 4, 3), xpd = NA)

image(
  x = seq_len(ncol(matriz_top5)),
  y = seq_len(nrow(matriz_top5)),
  z = t(matriz_top5[nrow(matriz_top5):1, , drop = FALSE]),
  col = colores_heatmap,
  zlim = c(-12, 12),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "Top 5 genes with the largest differences in quantification"
)

axis(2, at = 1:5, labels = rev(etiquetas_genes), las = 2, tick = FALSE, cex.axis = 0.75)

text(x = 1:3, y = rep(0.25, 3), labels = colnames(matriz_top5), srt = 45, adj = 1, cex = 0.85)

for (i in seq_len(nrow(matriz_top5))) {
  for (j in seq_len(ncol(matriz_top5))) {
    text(j, nrow(matriz_top5) - i + 1, labels = sprintf("%.2f", matriz_top5[i, j]), cex = 0.85)
  }
}

box()

mtext(
  paste(
    "Log2 ratio of Salmon estimated counts to featureCounts counts across three RNA-seq samples.",
    "Red (positive): higher Salmon quantification | Blue (negative): higher featureCounts quantification.",
    "Differences are consistent across samples; their underlying causes require further investigation.",
    sep = "\n"
  ),
  side = 1,
  line = 8.5,
  cex = 0.70,
  adj = 0.5
)

dev.off()


# 9. Anotación de genes discrepantes: Verificar la presencia de los cinco genes en ambas matrices y consultar la anotación Ensembl.
# Recuperar los nombres, biotipos y cromosomas utilizando el archivo GTF de GRCh38.115.

genes_interes <- top5_genes

verificacion_genes <- data.frame(Geneid = genes_interes, en_featureCounts = genes_interes %in% rownames(fc_counts), en_Salmon = genes_interes %in% rownames(salmon_counts))

print(verificacion_genes)

gtf_file <- "data/references/annotation/Homo_sapiens.GRCh38.115.gtf.gz"
lineas_genes <- system2("zgrep", args = c("-E", shQuote("^[^[:space:]]+[[:space:]]+[^[:space:]]+[[:space:]]+gene[[:space:]]"), shQuote(gtf_file)), stdout = TRUE)

campos_todos <- strsplit(lineas_genes, "\t", fixed = TRUE)

extraer_atributo <- function(texto, atributo) {
  patron <- paste0('.*', atributo, ' "([^"]+)".*')
  ifelse(grepl(paste0(atributo, ' "'), texto, fixed = TRUE), sub(patron, "\\1", texto), NA_character_)
}

genes_gtf <- data.frame(
  Geneid = sapply(campos_todos, function(x) extraer_atributo(x[9], "gene_id")),
  gene_name = sapply(campos_todos, function(x) extraer_atributo(x[9], "gene_name")),
  biotype = sapply(campos_todos, function(x) extraer_atributo(x[9], "gene_biotype")),
  chromosome = sapply(campos_todos, function(x) x[1]),
  start = as.numeric(sapply(campos_todos, function(x) x[4])),
  end = as.numeric(sapply(campos_todos, function(x) x[5]))
)

if (anyDuplicated(genes_gtf$Geneid)) stop("Se encontraron identificadores de genes duplicados en el GTF.")

indices_top5 <- match(top5_genes, genes_gtf$Geneid)

if (anyNA(indices_top5)) stop("No se encontraron todos los genes top 5 en la anotación GTF.")

anotacion_top5 <- genes_gtf[indices_top5, c("Geneid", "gene_name", "biotype", "chromosome"), drop = FALSE]

print(anotacion_top5)


# 10. Número de transcritos anotados: Contar los transcritos asociados a cada uno de los cinco genes discrepantes.
# Utilizar la tabla tx2gene derivada del mismo archivo FASTA de cDNA empleado para construir el índice de Salmon.

tx2gene <- readRDS("results/gene_expression/tx2gene_fasta.rds")

numero_transcritos <- sapply(genes_interes, function(gen) sum(tx2gene$GENEID == gen))

tabla_transcritos <- data.frame(Geneid = genes_interes, numero_transcritos = as.integer(numero_transcritos))

tabla_transcritos <- merge(tabla_transcritos, anotacion_top5, by = "Geneid", all.x = TRUE, sort = FALSE)
tabla_transcritos <- tabla_transcritos[match(top5_genes, tabla_transcritos$Geneid), ]

print(tabla_transcritos)


# 11. Visualización de transcritos anotados: Representar el número de transcritos asociados a los cinco genes discrepantes.
# Distinguir visualmente los genes protein-coding de los pseudogenes según la anotación GTF.

nombres_figura <- ifelse(is.na(tabla_transcritos$gene_name), "Processed pseudogene", tabla_transcritos$gene_name)
etiquetas_figura <- paste0(tabla_transcritos$Geneid, "\n(", nombres_figura, ")")

colores_figura <- ifelse(tabla_transcritos$biotype == "protein_coding", "#3B82B8", "#D98B45")

png("results/gene_expression/figures/top5_gene_transcript_counts.png", width = 1800, height = 1100, res = 200)

par(mar = c(5, 13, 5, 3), xpd = NA)

valores <- rev(tabla_transcritos$numero_transcritos)
etiquetas <- rev(etiquetas_figura)
colores <- rev(colores_figura)
posiciones <- seq_along(valores)

plot(
  valores,
  posiciones,
  xlim = c(0, max(29, valores + 5)),
  ylim = c(0.5, 5.5),
  yaxt = "n",
  pch = 19,
  cex = 2,
  col = colores,
  xlab = "Number of annotated transcripts",
  ylab = "",
  main = "Annotated transcripts in the five most discrepant genes",
  bty = "n"
)

axis(2, at = posiciones, labels = etiquetas, las = 1, tick = FALSE, cex.axis = 0.85)

segments(x0 = 0, y0 = posiciones, x1 = valores, y1 = posiciones, col = adjustcolor(colores, alpha.f = 0.35), lwd = 2)

points(valores, posiciones, pch = 19, cex = 2, col = colores)
text(valores + 0.8, posiciones, labels = valores, pos = 4, cex = 0.9)

legend("topright", inset = c(0, -0.12), legend = c("Protein-coding", "Pseudogene"), col = c("#3B82B8", "#D98B45"), pch = 19, pt.cex = 1.4, horiz = TRUE, bty = "n")

dev.off()


# 12. Evaluación del solapamiento genómico: Identificar genes anotados cuyas posiciones genómicas se solapan con los cinco genes discrepantes.
# Contar los genes solapados sin asumir que el solapamiento explica por sí solo las diferencias entre métodos.

posiciones_top5 <- genes_gtf[indices_top5, c("Geneid", "chromosome", "start", "end"), drop = FALSE]

print(posiciones_top5)

resultados_solapamiento <- do.call(rbind, lapply(seq_len(nrow(posiciones_top5)), function(i) {
  
  gen <- posiciones_top5[i, ]
  
  encontrados <- genes_gtf[
    genes_gtf$chromosome == gen$chromosome &
      genes_gtf$start <= gen$end &
      genes_gtf$end >= gen$start &
      genes_gtf$Geneid != gen$Geneid,
  ]
  
  data.frame(gen_principal = gen$Geneid, genes_solapados = nrow(encontrados))
}))

print(resultados_solapamiento)


# 13. Control de calidad de featureCounts: Revisar la asignación de lecturas a genes y las principales categorías de lecturas no asignadas.
# Calcular los porcentajes por muestra sin interpretar Unassigned_Ambiguity como equivalente a multimapping.

resumen_featurecounts <- read.delim("results/alignment_based/featurecounts/gene_counts_tail1.txt.summary", row.names = 1, check.names = FALSE)

print(resumen_featurecounts)

total_lecturas <- colSums(resumen_featurecounts)

porcentajes_featurecounts <- round(sweep(resumen_featurecounts, 2, total_lecturas, "/") * 100, 2)

print(porcentajes_featurecounts[
  c(
    "Assigned",
    "Unassigned_Unmapped",
    "Unassigned_NoFeatures",
    "Unassigned_Ambiguity",
    "Unassigned_MultiMapping"
  ),
])


# 14. Verificación de figuras: Comprobar que los tres archivos PNG se han guardado correctamente.
# Confirmar la existencia de las figuras de correlación, discrepancias y número de transcritos.

archivos_figuras <- file.path("results/gene_expression/figures", c("featurecounts_vs_salmon_3panels.png", "top5_gene_discrepancies_heatmap.png", "top5_gene_transcript_counts.png"))

data.frame(figura = basename(archivos_figuras), existe = file.exists(archivos_figuras))





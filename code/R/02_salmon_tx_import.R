# 1. Preparación del entorno: Cargar tximport y definir el directorio del proyecto para importar los resultados de Salmon.
# Este script utiliza las cuantificaciones de tres muestras y genera una matriz de expresión a nivel de gen.

library(tximport)
setwd("/Users/viancyverhoeven/OneDrive/javeriana/prueba_RNA-seq_analisis")

# 2. Identificación de muestras y archivos de Salmon: Definir las muestras y localizar sus archivos de cuantificación.
# Comprobar que los tres archivos quant.sf existen antes de continuar con la importación.

muestras <- c("SRR19020908", "SRR19020909", "SRR19020910")
archivos <- file.path("results/pseudoalignment/salmon_quant", muestras, "quant.sf")
names(archivos) <- muestras
stopifnot(all(file.exists(archivos)))

# 3. Identificación del transcriptoma de referencia: Localizar el archivo FASTA de Ensembl utilizado para construir el índice de Salmon.
# Utilizar los encabezados de este mismo archivo para mantener la correspondencia entre los identificadores de transcritos y genes.

archivo_fasta <- "data/references/transcriptome/Homo_sapiens.GRCh38.cdna.all.fa.gz"
stopifnot(file.exists(archivo_fasta))

# 4. Extracción de identificadores del FASTA: Recuperar los encabezados de las secuencias y extraer los identificadores de transcritos y genes.
# Conservar los identificadores completos de los transcritos, incluidas sus versiones, para que coincidan con Salmon.

encabezados <- system2("zgrep", c("'^>'", shQuote(archivo_fasta)), stdout = TRUE)
stopifnot(length(encabezados) > 0)
stopifnot(all(grepl(" gene:ENSG[0-9]+\\.[0-9]+", encabezados)))

ids_transcritos <- sub("^>(\\S+).*", "\\1", encabezados)
ids_genes <- sub(".* gene:(ENSG[0-9]+)\\.[0-9]+.*", "\\1", encabezados)
tx2gene_fasta <- data.frame(TXNAME = ids_transcritos, GENEID = ids_genes, stringsAsFactors = FALSE)

stopifnot(!anyDuplicated(tx2gene_fasta$TXNAME))
stopifnot(all(grepl("^ENSG[0-9]+$", tx2gene_fasta$GENEID)))
dim(tx2gene_fasta)

# 5. Validación de identificadores: Comparar los identificadores de transcritos de Salmon con los obtenidos del FASTA.
# Exigir una correspondencia completa en las tres muestras antes de realizar la agregación a nivel de gen.

correspondencia <- vapply(archivos, function(archivo) mean(read.delim(archivo)$Name %in% tx2gene_fasta$TXNAME) * 100, numeric(1))
print(correspondencia)
stopifnot(all(correspondencia == 100))

# 6. Importación y agregación con tximport: Importar las cuantificaciones de Salmon y resumir la expresión a nivel de gen.
# Los valores de counts son estimaciones de abundancia y pueden contener decimales, a diferencia de los conteos enteros de featureCounts.

txi <- tximport(archivos, type = "salmon", tx2gene = tx2gene_fasta)
dim(txi$counts)
colSums(txi$counts)

# 7. Validación de la matriz de expresión: Comprobar que la matriz contiene las tres muestras y valores válidos.
# Verificar que no existen valores ausentes, infinitos o negativos antes de guardar los resultados.

stopifnot(identical(colnames(txi$counts), muestras))
stopifnot(all(is.finite(txi$counts)))
stopifnot(all(txi$counts >= 0))

# 8. Exportación de resultados: Guardar el objeto tximport y la tabla de correspondencia entre transcritos y genes.
# Estos archivos permiten reutilizar los resultados en otros scripts sin repetir la importación de Salmon.

dir.create("results/gene_expression", recursive = TRUE, showWarnings = FALSE)
saveRDS(txi, "results/gene_expression/tximport_gene_counts.rds")
saveRDS(tx2gene_fasta, "results/gene_expression/tx2gene_fasta.rds")

stopifnot(file.exists("results/gene_expression/tximport_gene_counts.rds"))
stopifnot(file.exists("results/gene_expression/tx2gene_fasta.rds"))
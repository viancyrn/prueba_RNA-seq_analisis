
#!/usr/bin/env bash
set -euo pipefail

# FASE 2: CONVERSION DE SRA A FASTQ
# Dataset: GSE201955

# Ir al directorio principal del proyecto
cd "$(dirname "$0")/.."

mkdir -p 00_raw

while IFS= read -r SRR || [[ -n "$SRR" ]]; do

    # Ignorar líneas vacías y comentarios
    [[ -z "$SRR" || "$SRR" == \#* ]] && continue

    # Definir archivos de entrada y salida
    INPUT="00_raw/${SRR}.sra"
    OUTPUT="00_raw/${SRR}.fastq"

    # Comprobar que existe el archivo SRA
    if [[ ! -s "$INPUT" ]]; then
        echo "ERROR: No existe $INPUT"
        exit 1
    fi

    # Evitar convertir muestras ya procesadas
    if [[ -s "$OUTPUT" ]]; then
        echo "$SRR ya tiene FASTQ. Se omite."
        continue
    fi

    echo "Convirtiendo $SRR a FASTQ ..."

    # Convertir una muestra a la vez
    fasterq-dump "$INPUT" -O 00_raw -e 4 -p

    # Comprobar el resultado
    if [[ ! -s "$OUTPUT" ]]; then
        echo "ERROR: No se generó $OUTPUT"
        exit 1
    fi

    echo "$SRR convertido correctamente."

done < samples.txt

echo "Conversion finalizada."


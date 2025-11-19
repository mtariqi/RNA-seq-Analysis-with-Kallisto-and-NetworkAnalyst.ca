#!/bin/bash
# Student name: Md Tariqul Islam
# NUID: 002505253

# ---------------------- USER + COURSE DIRECTORIES ----------------------------------
# Location of course image

COURSE_DIR=$(readlink -f /courses/BINF6310.202610)
DATA_DIR="${COURSE_DIR}/data/RNA-Seq/rnaseq-mus-musculus-GSE240196"

# Kallisto Index file
INDEX="Mus_musculus.idx"

echo "Using directory: ${DIRECTORY}"
echo "Using index: ${INDEX}"

# ------------------ Loop through all GSM* files ---------------------------------


for FILE in $DATA_DIR/GSM*; do

    # Extract the sample name 
    BASENAME=$(basename "$FILE")

    # Create output directory per sample
    OUTPUT_DIR="kallisto-output-${BASENAME}"
    mkdir -p '${OUTPUT_DIR}'

    echo "================================"
    echo "processing sample: ${BASENAME}"
    echo "Input file: ${FILE}"
    echo "Output directory: ${OUTPUT_DIR}"



    # -------------- Run Kallisto via Apptainer ---------------------------------
    apptainer run \
        -B "/courses:/courses,/scratch:/scratch,/home:/home" \
        $COURSE_DIR/data/shared/binf6310_202530.sif \
        kallisto quant \
            -i ${INDEX} \
            -l 200 \
            -s 20 \
            -o "${OUTPUT_DIR}" \
            --single "${FILE}"

    echo "Finished ${BASENAME}"

done

echo "All samples processed!"

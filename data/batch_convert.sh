#!/bin/bash

# Check for exactly one argument: the target directory
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory_path>"
    echo "Example: $0 ../../clustering/modularity/DCHSBM/scenA1/"
    exit 1
fi

DIR="$1"

# Remove trailing slash if present
DIR="${DIR%/}"

# Loop over rep1 to rep25
for i in $(seq 1 25); do
    INPUT_TXT="${DIR}/rep${i}_he.txt"
    OUTPUT_HGR="${DIR}/rep${i}_he.hgr"
    GT_FILE="${DIR}/rep${i}_assign.txt"

    echo "Converting rep${i}..."
    python3 convert_to_hmetis.py "$INPUT_TXT" "$OUTPUT_HGR" "$GT_FILE"
done


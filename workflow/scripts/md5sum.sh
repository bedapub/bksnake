#!/bin/bash

# Directory containing the md5sum output files
CHECKSUM_FOLDER="$1"

# Check if folder is provided and exists
if [[ -z "$CHECKSUM_FOLDER" || ! -d "$CHECKSUM_FOLDER" ]]; then
    echo "Usage: $0 <folder_with_md5sum_files>"
    exit 1
fi

# Temporary associative array to store checksums and their occurrences
declare -A CHECKSUM_MAP

# Read all files in the folder
for FILE in "$CHECKSUM_FOLDER"/*; do
    # Skip if not a file
    [[ -f "$FILE" ]] || continue

    # Read each line of the file
    while read -r LINE; do
        # Extract the checksum (first field) from the line
        CHECKSUM=$(echo "$LINE" | awk '{print $}')
        
        # Skip empty or invalid lines
        [[ -n "$CHECKSUM" ]] || continue

        # Increment the count for this checksum
        ((CHECKSUM_MAP++))
    done < "$FILE"
done

# Check for duplicate checksums
DUPLICATES_FOUND=0
for CHECKSUM in "${!CHECKSUM_MAP[@]}"; do
    if [[ "${CHECKSUM_MAP[$CHECKSUM]}" -gt 1 ]]; then
        echo "Warning: Checksum '$CHECKSUM' appears in ${CHECKSUM_MAP[$CHECKSUM]} files!"
        DUPLICATES_FOUND=1
    fi
done

# Output final result
if [[ $DUPLICATES_FOUND -eq 0 ]]; then
    echo "All checksums are unique. No duplicates found."
else
    echo "Duplicate checksums detected. Please investigate the warnings above."
fi

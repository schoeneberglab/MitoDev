#!/bin/bash

# Usage: ./download.sh <destination_directory>
DEST_DIR="$1"

if [ -z "$DEST_DIR" ]; then
    echo "Usage: $0 <destination_directory>"
    exit 1
fi

# Create destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Hardcoded Google Drive link and filename
GDRIVE_LINK="https://drive.google.com/uc?id=1WNH7xLb9TAL5SojoPki8xnx9YS_EaS4D"
FILENAME="mitodev_data_sample.zip"

# Download using gdown
echo "Downloading $FILENAME..."
gdown "$GDRIVE_LINK" -O "${DEST_DIR}/${FILENAME}"

# Check if download succeeded
if [ ! -f "${DEST_DIR}/${FILENAME}" ]; then
    echo "Error: Download failed!"
    exit 1
fi

# Unzip the file
echo "Unzipping ${FILENAME}..."
unzip -o "${DEST_DIR}/${FILENAME}" -d "$DEST_DIR"

mkdir -p "./data/20231221 Gillian Lung Organoid/Sample 1/1/Processed_Data"
mv ./data/prod/* "./data/20231221 Gillian Lung Organoid/Sample 1/1/Processed_Data/"
rmdir ./data/prod

echo "Data extracted to $DEST_DIR"
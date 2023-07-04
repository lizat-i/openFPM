#!/bin/bash

# Create the new directory if it doesn't exist
dir="./modified_files"
if [[ ! -e $dir ]]; then
    mkdir $dir
elif [[ ! -d $dir ]]; then
    echo "$dir already exists but is not a directory" 1>&2
    exit 1
fi

# Iterate over all files in the current directory
for file in ./*; do
  # Check if the file name matches the pattern "Geometry_*"
  if [[ "${file##*/}" =~ Geometry_* ]]; then
    # Replace "Uint" with "UInt" in the file
    sed -i.bak 's/Uint/UInt/g' "$file"

    # Remove the backup file created by sed
    rm "${file}.bak"

    # Move the modified file to the new directory
    mv "$file" "$dir"
  fi
done

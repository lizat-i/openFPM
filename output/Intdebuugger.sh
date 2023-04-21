#!/bin/bash

# Iterate over all files in the current directory
for file in ./*; do
  # Check if the file name matches the pattern "myParticlesX"
  if [[ "${file##*/}" =~ ^Geometry.* ]]; then
    # Replace "Uint" with "UInt" in the file
    sed -i.bak 's/Uint/UInt/g' "$file"

    # Remove the backup file created by sed
    rm "${file}.bak"
  fi
done

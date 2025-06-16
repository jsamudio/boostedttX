#!/bin/bash

# Check if a filename is provided as an argument
if [ $# -eq 0 ]; then
  echo "Error: Please provide a filename as an argument."
  exit 1
fi

# Filename provided
file="$1"

# Use sed to modify the file in-place
sed -i 's/^/"file:/' "$file"  # Add "file:" at the beginning
sed -i 's/$/",/g' "$file"     # Add ", at the end (globally)

echo "Successfully modified '$file'"

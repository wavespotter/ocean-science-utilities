#!/bin/bash
OSU_VERSION="0.0.0"
# Run pdoc to generate documentation and redirect the output to README.pdf
doc_directory=docs/$OSU_VERSION
mkdir -p $doc_directory
mkdir -p $doc_directory/html
mkdir -p $doc_directory/md

pdoc src/ocean_science_utilities -o $doc_directory/html --html
pdoc src/ocean_science_utilities -o $doc_directory/md
pdoc src/ocean_science_utilities --pdf > $doc_directory/README.md

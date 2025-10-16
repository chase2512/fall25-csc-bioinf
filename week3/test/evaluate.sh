#!/bin/bash

# Print header
echo "Language    Runtime"
echo "-------------------"

# Run Python test
python3 test_phylo_python.py

# Run Codon test
codon run -release test_phylo_codon.py
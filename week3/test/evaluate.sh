#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="$(cd "${SCRIPT_DIR}/../code" && pwd)"

# Set PYTHONPATH to include the code directory
export PYTHONPATH="${CODE_DIR}:${PYTHONPATH}"

# Print header
echo "Language    Runtime"
echo "-------------------"

# Run Python test
python3 test_phylo_python.py

# Run Codon test (sys.path.insert works for Codon too)
codon run -release test_phylo_codon.py
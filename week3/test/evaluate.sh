#!/bin/bash
set -euo pipefail

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="$(cd "${SCRIPT_DIR}/../code" && pwd)"

# Set PYTHONPATH to include the code directory (for Python)
export PYTHONPATH="${CODE_DIR}:${PYTHONPATH:-}"

# Print header
echo "Language    Runtime"
echo "-------------------"

# Run Python test (uses biotite from pip)
python3 "${SCRIPT_DIR}/test_phylo_python.py"

# Run Codon test (uses our ported code)
# Codon doesn't support sys.path, so we create symlinks to the code files
cd "${SCRIPT_DIR}"

# Create symlinks if they don't exist
ln -sf ../code/tree_codon.py tree_codon.py 2>/dev/null || true
ln -sf ../code/upgma_codon.py upgma_codon.py 2>/dev/null || true
ln -sf ../code/nj_codon.py nj_codon.py 2>/dev/null || true

# Run the Codon test
codon run -release test_phylo_codon.py

# Clean up symlinks
rm -f tree_codon.py upgma_codon.py nj_codon.py
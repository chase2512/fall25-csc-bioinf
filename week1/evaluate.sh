#!/bin/bash
set -euo pipefail

# Function to calculate N50 from contig.fasta file
calculate_n50() {
    local fasta_file="$1"
    
    if [[ ! -f "$fasta_file" ]]; then
        echo "0"
        return
    fi
    
    # Extract sequence lengths and calculate N50
    python3 -c "
import sys

def calculate_n50(fasta_file):
    lengths = []
    with open(fasta_file, 'r') as f:
        current_length = 0
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_length > 0:
                    lengths.append(current_length)
                    current_length = 0
            else:
                current_length += len(line)
        if current_length > 0:
            lengths.append(current_length)
    
    if not lengths:
        return 0
    
    lengths.sort(reverse=True)
    total_length = sum(lengths)
    cumulative_length = 0
    
    for length in lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            return length
    
    return lengths[-1] if lengths else 0

print(calculate_n50('$fasta_file'))
"
}

# Function to format time in MM:SS format
format_time() {
    local seconds=$1
    printf "%d:%02d" $((seconds/60)) $((seconds%60))
}

# Use CODON_PYTHON from environment if available, otherwise auto-detect
if [[ -z "${CODON_PYTHON:-}" ]]; then
    # Try to find the correct Python library
    if [[ -f "/lib/x86_64-linux-gnu/libpython3.12.so.1.0" ]]; then
        export CODON_PYTHON="/lib/x86_64-linux-gnu/libpython3.12.so.1.0"
    elif [[ -f "/lib/x86_64-linux-gnu/libpython3.12.so" ]]; then
        export CODON_PYTHON="/lib/x86_64-linux-gnu/libpython3.12.so"
    else
        echo "Error: Could not find Python shared library"
        exit 1
    fi
fi

echo "Using CODON_PYTHON: ${CODON_PYTHON}"

# Verify the library exists
if [[ ! -f "${CODON_PYTHON}" ]]; then
    echo "Error: Python library not found at ${CODON_PYTHON}"
    echo "Available Python libraries:"
    find /lib /usr/lib -name "libpython*.so*" 2>/dev/null || echo "None found"
    exit 1
fi

# Change to the code directory
cd "$(dirname "$0")/code"

# Test datasets
datasets=("data1" "data2" "data3" "data4")

# Array to store results
declare -a results=()

for dataset in "${datasets[@]}"; do
    # Check if dataset exists in the data directory
    if [[ ! -d "../data/$dataset" ]]; then
        continue
    fi
    
    # Test Python version
    # Clean up any existing contig.fasta
    rm -f "../data/$dataset/contig.fasta"
    
    # Time the Python execution
    python_start=$(date +%s)
    if python3 main.py "../data/$dataset" >/dev/null 2>&1; then
        python_end=$(date +%s)
        python_runtime_seconds=$((python_end - python_start))
        python_runtime=$(format_time $python_runtime_seconds)
        python_n50=$(calculate_n50 "../data/$dataset/contig.fasta")
    else
        echo "Python execution failed for $dataset"
        python_runtime="FAILED"
        python_n50="N/A"
    fi
    
    # Store Python result
    results+=("$dataset python $python_runtime $python_n50")
    
    # Test Codon version
    # Clean up any existing contig.fasta
    rm -f "../data/$dataset/contig.fasta"
    
    # Time the Codon execution
    codon_start=$(date +%s)
    if codon run -release -plugin seq main_codon.py "../data/$dataset" >/dev/null 2>&1; then
        codon_end=$(date +%s)
        codon_runtime_seconds=$((codon_end - codon_start))
        codon_runtime=$(format_time $codon_runtime_seconds)
        codon_n50=$(calculate_n50 "../data/$dataset/contig.fasta")
    else
        echo "Codon execution failed for $dataset"
        codon_runtime="FAILED"
        codon_n50="N/A"
    fi
    
    # Store Codon result
    results+=("$dataset codon $codon_runtime $codon_n50")
done

# Print final table
echo "Dataset	Language	Runtime		N50"
echo "-------------------------------------------------------------------------------------------------------"

for result in "${results[@]}"; do
    read -r dataset language runtime n50 <<< "$result"
    printf "%-8s\t%-8s\t%-8s\t%s\n" "$dataset" "$language" "$runtime" "$n50"
done
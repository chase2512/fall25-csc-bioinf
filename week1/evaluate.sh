#!/bin/bash
set -euo pipefail

# Set unlimited stack size
ulimit -s unlimited

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

# Set up Codon PATH if not already available
if ! command -v codon &> /dev/null; then
    if [[ -f "${HOME}/.codon/bin/codon" ]]; then
        export PATH="${PATH}:${HOME}/.codon/bin"
        echo "Added Codon to PATH: ${HOME}/.codon/bin"
    fi
fi

# Use CODON_PYTHON from environment if available, otherwise auto-detect
if [[ -z "${CODON_PYTHON:-}" ]]; then
    if [[ -f "/lib/x86_64-linux-gnu/libpython3.12.so.1.0" ]]; then
        export CODON_PYTHON="/lib/x86_64-linux-gnu/libpython3.12.so.1.0"
    elif [[ -f "/lib/x86_64-linux-gnu/libpython3.12.so" ]]; then
        export CODON_PYTHON="/lib/x86_64-linux-gnu/libpython3.12.so"
    fi
fi

echo "Using CODON_PYTHON: ${CODON_PYTHON}"
echo "Current stack size: $(ulimit -s)"

# Change to the code directory
cd "$(dirname "$0")/code"

# Test only smaller datasets to avoid timeouts
datasets=("data1" "data2")

# Array to store results
declare -a results=()

for dataset in "${datasets[@]}"; do
    if [[ ! -d "../data/$dataset" ]]; then
        continue
    fi
    
    # Test Python version
    rm -f "../data/$dataset/contig.fasta"
    
    python_start=$(date +%s)
    if timeout 300 python3 main.py "../data/$dataset" >/dev/null 2>&1; then
        python_end=$(date +%s)
        python_runtime_seconds=$((python_end - python_start))
        python_runtime=$(format_time $python_runtime_seconds)
        python_n50=$(calculate_n50 "../data/$dataset/contig.fasta")
    else
        python_runtime="TIMEOUT"
        python_n50="N/A"
    fi
    
    results+=("$dataset python $python_runtime $python_n50")
    
    # Test Codon version
    rm -f "../data/$dataset/contig.fasta"
    
    codon_start=$(date +%s)
    if timeout 300 codon run -release -plugin seq main_codon.py "../data/$dataset" >/dev/null 2>&1; then
        codon_end=$(date +%s)
        codon_runtime_seconds=$((codon_end - codon_start))
        codon_runtime=$(format_time $codon_runtime_seconds)
        codon_n50=$(calculate_n50 "../data/$dataset/contig.fasta")
    else
        codon_runtime="TIMEOUT"
        codon_n50="N/A"
    fi
    
    results+=("$dataset codon $codon_runtime $codon_n50")
done

# Print final table
echo "Dataset	Language	Runtime		N50"
echo "-------------------------------------------------------------------------------------------------------"

for result in "${results[@]}"; do
    read -r dataset language runtime n50 <<< "$result"
    printf "%-8s\t%-8s\t%-8s\t%s\n" "$dataset" "$language" "$runtime" "$n50"
done
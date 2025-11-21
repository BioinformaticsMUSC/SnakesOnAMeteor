#!/bin/bash

# METEOR Snakemake Pipeline Test Runner
# This script sets up a test environment and runs the pipeline

set -e

echo "=== METEOR Snakemake Pipeline Test ==="

# Create test data directory structure
echo "Creating test environment..."
mkdir -p data/fastq
mkdir -p test_results

# Create minimal test FASTQ files (empty but with correct format)
echo "Creating test FASTQ files..."
cat > data/fastq/sample1_R1.fastq << 'EOF'
H4:
@seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@seq2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

cat > data/fastq/sample1_R2.fastq << 'EOF'
@seq1
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@seq2
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

# Compress the files
gzip -f data/fastq/*.fastq

# Create test configuration
echo "Creating test configuration..."
cat > config_test.yaml << 'EOF'
input_dir: "data/fastq"
output_dir: "test_results"
catalogue_name: "hs_10_4_gut"
cpus: 2
fast: true  # Use fast mode for testing
check_catalogue: false
minimum_memory_gb: 8
EOF

echo "Test environment created!"
echo ""
echo "To run the pipeline:"
echo "  snakemake --configfile config_test.yaml --cores 2 --use-conda --dry-run"
echo ""
echo "To run for real:"
echo "  snakemake --configfile config_test.yaml --cores 2 --use-conda"
echo ""
echo "Files created:"
echo "  - data/fastq/sample1_R1.fastq.gz"
echo "  - data/fastq/sample1_R2.fastq.gz"
echo "  - config_test.yaml"
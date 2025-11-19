# METEOR Shotgun Metagenomics Snakemake Pipeline

A Snakemake pipeline for quantitative metagenomics profiling using METEOR, converted from the original Nextflow pipeline.

## Overview

METEOR is a platform for quantitative metagenomics profiling of complex ecosystems. It performs:
- Species-level taxonomic profiling (Bacteria, Archaea, Eukaryotes)
- Functional analysis
- Strain-level population structure inference

## Requirements

- Snakemake (≥6.0)
- Conda/Mamba
- Python ≥3.11

## Installation

1. Clone or download this pipeline
2. Install Snakemake:
   ```bash
   conda install -c conda-forge -c bioconda snakemake
   ```

## Usage

### 1. Prepare your data

Place paired-end FASTQ files in a directory with the naming pattern:
- `*_R1*.fastq.gz` and `*_R2*.fastq.gz`
- Or `*_R1*.fastq` and `*_R2*.fastq`

### 2. Configure the pipeline

Edit `config.yaml`:
- Set `input_dir` to your FASTQ directory
- Set `output_dir` for results
- Choose a `catalogue_name` or provide a custom `catalogue_path`
- Adjust `cpus`, `fast` mode, and memory settings as needed

### 3. Run the pipeline

```bash
# Dry run to check the workflow
snakemake --configfile config.yaml --dry-run

# Run with 8 cores
snakemake --configfile config.yaml --cores 8 --use-conda

# Run on a cluster (example with SLURM)
snakemake --configfile config.yaml --cores 8 --use-conda --cluster "sbatch -t 60 -c {threads}"
```

## Pipeline Steps

1. **Download/Prepare Catalogue**: Downloads prebuilt microbial gene catalogue
2. **FASTQ Processing**: Indexes and prepares FASTQ files
3. **Mapping**: Maps reads against the gene catalogue
4. **Profiling**: Generates taxonomic and functional profiles
5. **Merging**: Combines profiles from all samples
6. **Strain Analysis**: Performs strain-level analysis (if not in fast mode)
7. **Tree Building**: Constructs phylogenetic trees from strain data

## Output

Results are written to the specified `output_dir`:
- `merged/`: Combined taxonomic and functional profiles
- `tree/`: Phylogenetic trees (if strain analysis performed)
- `meteor_report.html`: Summary report

## Configuration Options

### Catalogue Options
Choose from prebuilt catalogues:
- `hs_10_4_gut`: Human gut microbiome
- `hs_8_4_oral`: Human oral microbiome  
- `hs_2_9_skin`: Human skin microbiome
- `fc_1_3_gut`: Cat gut microbiome
- `clf_1_0_gut`: Dog gut microbiome
- `mm_5_0_gut`: Mouse gut microbiome
- `ssc_9_3_gut`: Pig gut microbiome
- And others (see config.yaml)

### Processing Modes
- **Normal mode**: Full taxonomic + functional analysis
- **Fast mode** (`fast: true`): Taxonomic analysis only, reduced memory usage

## Memory Requirements

- Normal mode: ~30GB+ RAM (depends on sample size)
- Fast mode: ~10GB RAM maximum
- Memory scales with number of reads per sample

## Troubleshooting

### Common Issues

1. **No FASTQ files found**: Check file naming pattern and input directory
2. **Memory errors**: Increase `minimum_memory_gb` or enable fast mode
3. **Catalogue download fails**: Check internet connection and catalogue name

### Cleaning Up

Remove intermediate files:
```bash
snakemake clean --configfile config.yaml
```

## Comparison with Original Nextflow Pipeline

This Snakemake version maintains the same functionality as `nf-meteor.nf`:

| Feature | Nextflow | Snakemake |
|---------|----------|-----------|
| Catalogue download | ✓ | ✓ |
| FASTQ processing | ✓ | ✓ |
| Read mapping | ✓ | ✓ |
| Taxonomic profiling | ✓ | ✓ |
| Functional profiling | ✓ | ✓ |
| Profile merging | ✓ | ✓ |
| Strain analysis | ✓ | ✓ |
| Tree construction | ✓ | ✓ |
| Memory management | ✓ | ✓ |
| Fast mode | ✓ | ✓ |

## Citation

If you use this pipeline, please cite:
- METEOR: Ghozlane et al. "Accurate profiling of microbial communities for shotgun metagenomic sequencing with Meteor2." Microbiome (2025)
- Snakemake: Köster & Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics (2012)
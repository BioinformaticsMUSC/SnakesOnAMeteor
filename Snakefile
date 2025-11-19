#!/usr/bin/env python3
"""
METEOR Shotgun Metagenomics Snakemake Pipeline

A Snakemake pipeline for quantitative metagenomics profiling using METEOR.
Converted from the original Nextflow pipeline (nf-meteor.nf).

Usage:
    snakemake --configfile config.yaml --cores 8
"""

import os
import glob
from pathlib import Path

# Load configuration
configfile: "config.yaml"

# Validate required parameters
if not config.get("input_dir"):
    raise ValueError("ERROR: input_dir must be specified in config.yaml")
if not config.get("output_dir"):
    raise ValueError("ERROR: output_dir must be specified in config.yaml")
if not config.get("catalogue_name") and not config.get("catalogue_path"):
    raise ValueError("ERROR: Either catalogue_name or catalogue_path must be specified in config.yaml")

# Set default values
INPUT_DIR = config["input_dir"]
OUTPUT_DIR = config["output_dir"]
CATALOGUE_NAME = config.get("catalogue_name", "")
CATALOGUE_PATH = config.get("catalogue_path", "")
CPUS = config.get("cpus", 4)
FAST_MODE = config.get("fast", False)
CHECK_CATALOGUE = config.get("check_catalogue", False)
MINIMUM_MEMORY_GB = config.get("minimum_memory_gb", 30)

# Allowed catalogue names
ALLOWED_CATALOGUES = [
    "fc_1_3_gut", "gg_13_6_caecal", "clf_1_0_gut", "hs_10_4_gut", 
    "hs_8_4_oral", "hs_2_9_skin", "mm_5_0_gut", "oc_5_7_gut", 
    "rn_5_9_gut", "ssc_9_3_gut"
]

# Validate catalogue name
if CATALOGUE_NAME and CATALOGUE_NAME not in ALLOWED_CATALOGUES:
    raise ValueError(f"ERROR: Invalid catalogue name '{CATALOGUE_NAME}'. "
                    f"Allowed catalogues are: {', '.join(ALLOWED_CATALOGUES)}")

# Find paired-end FASTQ files
def get_samples():
    """Find paired-end FASTQ files in the input directory."""
    patterns = [
        f"{INPUT_DIR}/*_R1*.fastq.gz",
        f"{INPUT_DIR}/*_R1*.fastq",
        f"{INPUT_DIR}/*_R1*.fq.gz",
        f"{INPUT_DIR}/*_R1*.fq"
    ]
    
    samples = {}
    for pattern in patterns:
        for r1_file in glob.glob(pattern):
            # Extract sample name
            basename = os.path.basename(r1_file)
            sample_name = basename.replace("_R1", "").split(".")[0]
            
            # Find corresponding R2 file
            r2_file = r1_file.replace("_R1", "_R2")
            if os.path.exists(r2_file):
                samples[sample_name] = {"R1": r1_file, "R2": r2_file}
    
    if not samples:
        raise ValueError(f"No paired-end FASTQ files found in {INPUT_DIR}")
    
    return samples

SAMPLES = get_samples()

# Determine catalogue directory based on configuration
if CATALOGUE_NAME:
    CATALOGUE_DIR = f"{CATALOGUE_NAME}{'_taxo' if FAST_MODE else ''}"
else:
    CATALOGUE_DIR = CATALOGUE_PATH

# Rule priorities
ruleorder: meteor_download > meteor_custom_catalogue

# Target rule
rule all:
    input:
        f"{OUTPUT_DIR}/merged",
        f"{OUTPUT_DIR}/tree" if not FAST_MODE else []

# Download prebuilt catalogue
rule meteor_download:
    output:
        directory(f"{CATALOGUE_NAME}{'_taxo' if FAST_MODE else ''}")
    params:
        catalogue_name = CATALOGUE_NAME,
        fast_option = "--fast" if FAST_MODE else "",
        check_option = "-c" if CHECK_CATALOGUE else ""
    conda:
        "envs/meteor.yaml"
    shell:
        """
        meteor download -i {params.catalogue_name} -o ./ {params.check_option} {params.fast_option}
        """

# Use custom catalogue (alternative to download)
rule meteor_custom_catalogue:
    input:
        CATALOGUE_PATH if CATALOGUE_PATH else []
    output:
        directory("custom_catalogue")
    run:
        if CATALOGUE_PATH:
            import shutil
            shutil.copytree(input[0], output[0])

# Process FASTQ files
rule meteor_fastq:
    input:
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        directory("results/{sample}/fastq")
    params:
        sample_dir = "results/{sample}/input"
    conda:
        "envs/meteor.yaml"
    shell:
        """
        # Create input directory and copy FASTQ files
        mkdir -p {params.sample_dir}
        cp {input.r1} {params.sample_dir}/
        cp {input.r2} {params.sample_dir}/
        
        # Run meteor fastq
        meteor fastq -i {params.sample_dir} -p -o {output}
        """

# Mapping reads to catalogue
rule meteor_mapping:
    input:
        fastq_dir = "results/{sample}/fastq",
        catalogue = CATALOGUE_DIR
    output:
        directory("results/{sample}/mapping")
    params:
        cpus = CPUS
    resources:
        mem_mb = lambda wildcards: max(MINIMUM_MEMORY_GB * 1024, 10 * 1024 if FAST_MODE else MINIMUM_MEMORY_GB * 1024)
    threads: CPUS
    conda:
        "envs/meteor.yaml"
    shell:
        """
        meteor mapping -i {input.fastq_dir} -r {input.catalogue} -t {threads} -o {output} --kf
        """

# Generate taxonomic and functional profiles
rule meteor_profile:
    input:
        mapping_dir = "results/{sample}/mapping",
        catalogue = CATALOGUE_DIR
    output:
        directory("results/{sample}/profile")
    threads: CPUS
    conda:
        "envs/meteor.yaml"
    shell:
        """
        meteor profile -i {input.mapping_dir} -r {input.catalogue} -o {output}
        """

# Merge profiles from all samples
rule meteor_merge:
    input:
        profiles = expand("results/{sample}/profile", sample=SAMPLES.keys()),
        catalogue = CATALOGUE_DIR
    output:
        directory(f"{OUTPUT_DIR}/merged")
    params:
        profile_parent_dir = "results"
    resources:
        mem_mb = lambda wildcards, input: max(2048, len(input.profiles) * 200 * 1.5)
    conda:
        "envs/meteor.yaml"
    shell:
        """
        # Create temporary directory with all profiles
        mkdir -p temp_profiles
        for profile_dir in {input.profiles}; do
            cp -r "$profile_dir" temp_profiles/
        done
        
        meteor merge -i temp_profiles -r {input.catalogue} -o {output} -s
        
        # Clean up
        rm -rf temp_profiles
        """

# Strain profiling (only in non-fast mode)
rule meteor_strain:
    input:
        mapping_dir = "results/{sample}/mapping",
        catalogue = CATALOGUE_DIR
    output:
        directory("results/{sample}/strain")
    resources:
        mem_mb = lambda wildcards: max(MINIMUM_MEMORY_GB * 1024, 10 * 1024 if FAST_MODE else MINIMUM_MEMORY_GB * 1024)
    conda:
        "envs/meteor.yaml"
    shell:
        """
        meteor strain -i {input.mapping_dir} -r {input.catalogue} -o {output}
        """

# Build phylogenetic tree from strains
rule meteor_tree:
    input:
        strains = expand("results/{sample}/strain", sample=SAMPLES.keys())
    output:
        directory(f"{OUTPUT_DIR}/tree")
    threads: CPUS
    conda:
        "envs/meteor.yaml"
    shell:
        """
        # Create temporary directory with all strains
        mkdir -p temp_strains
        for strain_dir in {input.strains}; do
            if [ -d "$strain_dir" ]; then
                cp -r "$strain_dir"/* temp_strains/ 2>/dev/null || true
            fi
        done
        
        # Only run tree if we have strain data
        if [ "$(ls -A temp_strains)" ]; then
            meteor tree -i temp_strains -r -o {output} -t {threads}
        else
            mkdir -p {output}
            echo "No strain data available for tree construction" > {output}/no_strains.txt
        fi
        
        # Clean up
        rm -rf temp_strains
        """

# Clean up intermediate files
rule clean:
    shell:
        """
        rm -rf results/*/fastq results/*/mapping results/*/profile results/*/strain
        """

# Generate summary report
rule report:
    input:
        merged = f"{OUTPUT_DIR}/merged",
        tree = f"{OUTPUT_DIR}/tree"
    output:
        f"{OUTPUT_DIR}/meteor_report.html"
    run:
        import datetime
        
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>METEOR Pipeline Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .section {{ margin: 20px 0; }}
                .file-list {{ background-color: #f9f9f9; padding: 15px; border-radius: 5px; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>METEOR Metagenomics Pipeline Report</h1>
                <p>Generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p>Catalogue: {CATALOGUE_NAME or CATALOGUE_PATH}</p>
                <p>Fast mode: {FAST_MODE}</p>
                <p>Samples processed: {len(SAMPLES)}</p>
            </div>
            
            <div class="section">
                <h2>Sample Information</h2>
                <div class="file-list">
                    <ul>
                    {"".join(f"<li>{sample}: {info['R1']} + {info['R2']}</li>" for sample, info in SAMPLES.items())}
                    </ul>
                </div>
            </div>
            
            <div class="section">
                <h2>Output Files</h2>
                <div class="file-list">
                    <p><strong>Merged profiles:</strong> {OUTPUT_DIR}/merged/</p>
                    {"<p><strong>Phylogenetic tree:</strong> " + OUTPUT_DIR + "/tree/</p>" if not FAST_MODE else ""}
                </div>
            </div>
            
            <div class="section">
                <h2>Pipeline Configuration</h2>
                <div class="file-list">
                    <p>CPUs used: {CPUS}</p>
                    <p>Minimum memory: {MINIMUM_MEMORY_GB} GB</p>
                    <p>Check catalogue: {CHECK_CATALOGUE}</p>
                </div>
            </div>
        </body>
        </html>
        """
        
        with open(output[0], 'w') as f:
            f.write(html_content)
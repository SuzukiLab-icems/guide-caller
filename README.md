README for guide-caller v1.0.0
Overview
guide-caller v1.0.0 is a comprehensive toolset designed for processing and analyzing GeCKO v2 CRISPR Knockout Pooled Library screening. This package includes two primary components: guide_caller.sh, a shell script that orchestrates the data processing workflow, and matrix_shaper.py, a Python script for advanced data analysis and matrix generation.
Installation
To use guide-caller, ensure you have the following prerequisites installed:
A UNIX-like environment (for running shell scripts)
Python 3.9.16 (developer's environment)
Required Python libraries: NumPy, pandas
Cutadapt v4.1 (developer's environment)
FastQC v0.12.1 (developer's environment)
MAGeCK v0.5.9.4 (developer's environment)
After ensuring the prerequisites, download the guide-caller package and extract its contents to your desired directory.
Usage
guide_caller.sh
This script processes FASTQ files and extracts and counts sgRNA sequences. It uses various genomic data processing tools and integrates with matrix_shaper.py.

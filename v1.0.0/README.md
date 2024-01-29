# **guide-caller v1.0.0**
## Overview
guide-caller v1.0.0 is a comprehensive toolset designed for processing and analyzing GeCKO v2 CRISPR Knockout Pooled Library screening.
This package includes two primary components: guide_caller.sh, a shell script that orchestrates the data processing workflow, 
and matrix_shaper.py, a Python script for advanced data analysis and matrix generation.

## Installation
To use guide-caller, ensure you have the following prerequisites installed:
•	A UNIX-like environment (for running shell scripts)
•	Python 3.9.16 (developer's environment)
•	Required Python libraries: NumPy, pandas
•	Cutadapt v4.1 (developer's environment)
•	FastQC v0.12.1 (developer's environment)
•	MAGeCK v0.5.9.4 (developer's environment)
After ensuring the prerequisites, download the guide-caller package and extract its contents to your desired directory.

## Usage
### guide_caller.sh
This script processes FASTQ files and extracts and counts sgRNA sequences. It uses various genomic data processing tools and integrates with matrix_shaper.py.

### usage:
sh guide_caller.sh [OPTIONS]
Options:
  -h          Display help
  -i VALUE    Specify the input directory (ex. "HN00171248" etc...)
  -f VALUE    Specify the alignment file (ex. Human_GeCKOv2_Library.csv)
  -c VALUE    Specify the CPU core

### matrix_shaper.py
A Python script for transforming and analyzing count data from the CRISPR experiments. It generates a result matrix with statistical analyses.

### usage:
python matrix_shaper.py -i <count_file> -o <output_dir>

## Authors
Noguchi Yuki (Jun Suzuki lab)
Contact: nyuhki21@gmail.com, jsuzuki@icems.kyoto-u.ac.jp

## Citation
If you use this tool in your research, please cite:
Noguchi, Y., Onodera, Y., Maruoka, M., Miyamoto, T., Kosako, H., Suzuki., J. 2024. In vivo CRISPR screening directly targeting testicular cells. Cell Genomics.

## License
This software is released under the MIT License. See LICENSE for more details.

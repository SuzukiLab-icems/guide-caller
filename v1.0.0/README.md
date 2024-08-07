# guide-caller v1.0.0
## Overview
`guide-caller v1.0.0` is a comprehensive toolset designed for processing and analyzing GeCKO v2 CRISPR Knockout Pooled Library screening.
This package includes two primary components: guide_caller.sh, a shell script that orchestrates the data processing workflow, 
and matrix_shaper.py, a Python script for advanced data analysis and matrix generation.

## Requirements
To use guide-caller, ensure you have the following prerequisites installed according to `mageck.yaml`.
After ensuring the prerequisites, download the guide-caller package and extract its contents to your desired directory.

## Usage
### guide_caller.sh
This script processes FASTQ files and extracts and counts sgRNA sequences. It uses various genomic data processing tools and integrates with `matrix_shaper.py`.
Before run this script, please check your directory architechture by `tree guide-caller`

```
guide-caller
|
|-v1.0.0
|	|-guide_caller.sh
|	|-matrix_shaper.py
|
|-alignment_files <- *Do not change the name!
|	|-alignment_file.csv
|
|-your directory
	|-sample_name
		|-sample_name.fastq <- fastq format should be used (not fq, *,gz)
	|-•••
		|-•••
```

### usage:
```bash
cd ./guide-caller

sh ./v1.0.0/guide_caller.sh [OPTIONS]
Options:
  -h          Display help
  -i VALUE    Specify the input directory (ex. "HN00171248" etc...)
  -f VALUE    Specify the alignment file (ex. Human_GeCKOv2_Library.csv)
  -c VALUE    Specify the CPU core

ex): 'sh ./v1.0.0/guide_caller.sh -i <your directory> -f <alignment_file.csv> -c <CPU core>'
```

### matrix_shaper.py
A Python script for transforming and analyzing count data from the CRISPR experiments. It generates a result matrix with statistical analyses in `guide_caller.sh`.

### usage:
```bash
python matrix_shaper.py -i <count_file> -o <output_dir>
```

## Authors
Noguchi Yuki (Jun Suzuki lab)
Contact: nyuhki21@gmail.com, jsuzuki@icems.kyoto-u.ac.jp

## Citation
If you use this tool in your research, please cite:
Noguchi, Y., Onodera, Y., Maruoka, M., Miyamoto, T., Kosako, H., Suzuki., J. 2024. In vivo CRISPR screening directly targeting testicular cells. Cell Genomics.

## License
This software is released under the MIT License. See LICENSE for more details.

## References
1. Li, W., Xu, H., Xiao, T., Cong, L., Love, M.I., Zhang, F., Irizarry, R.A., Liu, J.S., Brown, M.A., and Liu, X.S. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology. 2014; 15: 554. 10.1186/s13059-014-0554-4
2. Martin, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal. 2011; 17: 10-12. 10.14806/ej.17.1.200
3. FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

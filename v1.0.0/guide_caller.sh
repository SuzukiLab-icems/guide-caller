#!/usr/sh

###########################################################################
#guide_caller_v1.0.0/guide_caller.sh
#
#	 Copyright (c) 2024, Noguchi Yuki (Jun Suzuki lab)
#	 This software is released under the MIT License, see LICENSE (https://opensource.org/license/mit/).
#    @citation: Noguchi, Y., Onodera, Y., Maruoka, M., Miyamoto, T., Kosako, H., Suzuki., J. 2024. In vivo CRISPR screening directly targeting testicular cells. Cell Genomics.
#    @author:  Noguchi Yuki
#    @contact: nyuhki21@gmail.com,jsuzuki@icems.kyoto-u.ac.jp
#
##REFERENCE:
#1.	Li, W., Xu, H., Xiao, T., Cong, L., Love, M.I., Zhang, F., Irizarry, R.A., Liu, J.S., Brown, M.A., and Liu, X.S. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology. 2014; 15: 554. 10.1186/s13059-014-0554-4
#2. Martin, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal. 2011; 17: 10-12. 10.14806/ej.17.1.200
#3. FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#
##Detailed Usage:
#You can run guide_caller by typing 'sh guide_caller_v1.0.0/guide_caller.sh -i <your directory> -f <alignment_file.csv> -c <CPU core>' in 'guide_caller' directory. Please type pwd and confirm your current directory is 'guide_caller.'
#
#[Directory Architecture]
#guide_caller
#	|
#	|-v1.0.0
#	|	|-guide_caller.sh
#	|	|-matrix_shaper.py
#	|
#	|-alignment_files <- *Do not change the name!
#	|	|-alignment_file.csv
#	|
#	|-your directory
#		|-sample_name
#			|-sample_name.fastq <- *.fastq should be used!
#		|-•••
#			|-•••
#If you want to change the directory, pls fix line.68 according to your preference.
###########################################################################

input_directory=""
alignment_file=""
cpu=""

function usage {
    cat <<EOM
Usage: $(basename "$0") [OPTION]...
    -h          Display help
    -i VALUE    Specify the input directory (ex. "HN00171248" etc...)
    -f VALUE    Specify the alignment_file (ex.Human_GeCKOv2_Library.csv)
    -c VALUE	Specify the CPU core
EOM
    exit 2
}

function argparse {
    while getopts ":i:f:c:h" opt; do
        case "$opt" in
			i) input_directory="${OPTARG}" ;;
            f) alignment_file="${OPTARG}" ;;
			c) cpu="${OPTARG}" ;;
            h) usage ;;
            \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
            :)  echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
        esac
    done
    shift "$((OPTIND - 1))"
}

function config {
	dir=$(pwd)"/v1.0.0"
	echo "###########################################################"
	echo "Start guide_caller.."
	echo " "
	echo "///Configuration///"
	echo "version==========:guide_caller_v1.0.0"
	echo "package directory:${dir}"
	echo "input directory==:${input_directory}"
	echo "alignment file===:${alignment_file}"
	echo "CPU core ========:${cpu}"
	echo " "
	echo "###########################################################"
}

function env_setting {
	LIST=$(ls ./${input_directory})
	for SampleID in `echo ${LIST}`;do
		mkdir ./${input_directory}/${SampleID}/${SampleID}_fastqc
		mkdir ./${input_directory}/${SampleID}/mageck_result
	done
}

function sgRNA_extraction {
    FILE=$(cat ./${input_directory}/${SampleID}/${SampleID}.fastq \
        | head -n 2 \
        | tail -n 1 \
        | cut -c 1-5
    )
    if [ "${FILE}" = "ATCGA" ]; then
        cutadapt \
            -j ${cpu} \
            -u 30 \
            ./${input_directory}/${SampleID}/${SampleID}.fastq > ./${input_directory}/${SampleID}/5-trimmed_${SampleID}.fastq
        cutadapt \
            -j ${cpu} \
            -u -1 \
			./${input_directory}/${SampleID}/5-trimmed_${SampleID}.fastq > ./${input_directory}/${SampleID}/trimmed_${SampleID}.fastq
    elif [ "${FILE}" = "TAGCT" ]; then
        cutadapt \
            -j ${cpu} \
            -u 29 \
            ./${input_directory}/${SampleID}/${SampleID}.fastq > ./${input_directory}/${SampleID}/5-trimmed_${SampleID}.fastq
        cutadapt \
            -j ${cpu} \
            -u -2 \
			./${input_directory}/${SampleID}/5-trimmed_${SampleID}.fastq > ./${input_directory}/${SampleID}/trimmed_${SampleID}.fastq
    elif [ "${FILE}" = "CAAGT" ]; then
        cutadapt \
            -j ${cpu} \
            -u 28 \
            ./${input_directory}/${SampleID}/${SampleID}.fastq > ./${input_directory}/${SampleID}/5-trimmed_${SampleID}.fastq
        cutadapt \
            -j ${cpu} \
            -u -3 \
			./${input_directory}/${SampleID}/5-trimmed_${SampleID}.fastq > ./${input_directory}/${SampleID}/trimmed_${SampleID}.fastq
    else
        cutadapt \
            -j ${cpu} \
            -u 31 \
            ./${input_directory}/${SampleID}/${SampleID}.fastq > ./${input_directory}/${SampleID}/trimmed_${SampleID}.fastq
	fi
	fastqc --nogroup -o ./${input_directory}/${SampleID}/${SampleID}_fastqc ./${input_directory}/${SampleID}/trimmed_${SampleID}.fastq
}

function sgRNA_count {
    mageck count \
        -l ./alignment_files/${alignment_file} \
        -n ${SampleID} \
        --fastq ./${input_directory}/${SampleID}/trimmed_${SampleID}.fastq
	mv ${SampleID}.* ./${input_directory}/${SampleID}/mageck_result
	mv *.Rnw *.R ./${input_directory}/${SampleID}/mageck_result
}

function main {
    argparse "$@"
    config
	env_setting
	for SampleID in `echo ${LIST}`;do
		sgRNA_extraction
		sgRNA_count
		python ${dir}/matrix_shaper.py -o ${input_directory}/${SampleID}/mageck_result/${SampleID} \
		< ${input_directory}/${SampleID}/mageck_result/${SampleID}.count.txt
	done
}

main "$@"
rm ${input_directory}/*/trimmed_*.fastq ${input_directory}/*/5-trimmed_*.fastq
echo "done."
exit=0

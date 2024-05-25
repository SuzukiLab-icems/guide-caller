#!/usr/zsh
echo "guide-caller main program"
echo "Copyright (c) 2023 Yuki NOGUCHI, Jun Suzuki lab"
echo "This code is free software; you can redistribute it and/or modify it"
echo "under the terms of the BSD License (see the file COPYING included with"
echo "the distribution)."
echo "@version: v3.0.0"
echo "@author:  Yuki NOGUCHI"
echo "@contact: nyuhki21@gmail.com / jsuzuki@icems.kyoto-u.ac.jp"
echo "CMD: sh guide-caller/v3.0.0/main.sh -m <matrix_table>"

meta_table=""
VERSION="v3.0.0"
FILE_PATH=./${pwd}/${VERSION}

function usage {
    cat <<EOM
Usage: $(basename "$0") [OPTION]...
    -h          Display help
    -m VALUE	Specify the meta table (ex. sample_info.csv)
EOM
    exit 2
}

function argparse {
	while getopts ":m:h" opt; do
		case "$opt" in
			m) meta_table="${OPTARG}" ;;
			h) usage ;;
			\?) echo "Invalid option: -$OPTARG" >&2; usage ;;
			:)  echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
		esac
	done
	shift "$((OPTIND - 1))"
}

function config {
	line=$(echo "$line" | tr -d '\r')
	library=$(awk -F ',' '{print $9}' <<< "$line")
	fastq_directory=$(awk -F ',' '{print $2}' <<< "$line")
	sample=$(awk -F ',' '{print $3}' <<< "$line")
	alignment=$(awk -F ',' '{print $4}' <<< "$line")
	alignment=".${pwd}/alignment_files/${alignment/#\~/$HOME}"
	barcode_1=$(awk -F ',' '{print $5}' <<< "$line")
	barcode_2=$(awk -F ',' '{print $6}' <<< "$line")
	barcode_length=$(awk -F ',' '{print $7}' <<< "$line")
	cpu_core=$(awk -F ',' '{print $8}' <<< "$line")
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	echo "Start analyzing ${sample}.."
	echo "Library type: ${library}"
	echo "Directory of fastq file: ${fastq_directory}"
	echo "Path of alignment file: ${alignment}"
	echo "barcode_fwd:${barcode_1}"
	echo "barcode_rev:${barcode_2}"
	echo "barcode_length:${barcode_length}"
	echo "CPU:${cpu_core}"
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
}

function exec_CRISPRko {
	sh ${FILE_PATH}/src/trimCRISPR-ko.sh ${fastq_directory} ${sample} ${barcode_length} ${cpu_core}
	sh ${FILE_PATH}/src/countGuide-ko.sh ${fastq_directory} ${sample} ${alignment} ${FILE_PATH}
}

function exec_CRISPRia {
	sh ${FILE_PATH}/src/trimCRISPR-ia.sh ${fastq_directory} ${sample} ${barcode_length} ${cpu_core}
	sh ${FILE_PATH}/src/countGuide-ia.sh ${fastq_directory} ${sample} ${alignment} ${FILE_PATH}
}

function main {
	argparse "$@"
	tail -n+2 ${meta_table} | while read line; do
		config
		if [ "${library}" = 'ko' ];then
			exec_CRISPRko "$@"
		elif [ "${library}" = 'ia' ];then
			exec_CRISPRia "$@"
		fi
	done
	rm -rf ./*/*/*tmp.fastq.gz ./*/*/*trimmed_1.fastq.gz ./*/*/*trimmed_2.fastq.gz
}

echo "Start guide caller version ${VERSION}.."
main "$@"
echo "Done."


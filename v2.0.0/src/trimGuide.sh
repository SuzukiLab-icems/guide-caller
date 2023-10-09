#!/usr/sh
echo "guide-caller trimGuide program"
echo "Copyright (c) 2023 Yuki NOGUCHI, Jun Suzuki lab"
echo "This code is free software; you can redistribute it and/or modify it"
echo "under the terms of the BSD License (see the file COPYING included with"
echo "the distribution)."
echo "@version: v2.0.0"
echo "@author:  Yuki NOGUCHI"
echo "@contact: jsuzuki AT icems.kyoto-u.ac.jp"

DIR=$1
i=$2
b=$3
c=$4
VERSION="v2.0.0"
FILE_PATH=~/package/guide-caller/${VERSION}

echo "~~~~~~~~~~~~Welcome to trimGuide.sh~~~~~~~~~~~~~"

tail -n+2 ${i} | while read line; do
	sample=$(awk -F ',' '{print $2}' <<< "$line")
	alignment=$(awk -F ',' '{print $3}' <<< "$line")
	direction=$(awk -F ',' '{print $4}' <<< "$line")
	barcode_1=$(awk -F ',' '{print $5}' <<< "$line")
	barcode_2=$(awk -F ',' '{print $6}' <<< "$line")
	length=$(awk -F ',' '{print $7}' <<< "$line")
	if [ "${b}" = "barcode" ]; then
		trimming="barcode"
	else
		trimming="default"
	fi
	echo "Start analyzing your sample"
	echo "sample:${sample}"
	echo "alignment file:${alignment}"
	echo "trim mode:${trimming}"
	echo "barcode_fwd:${barcode_1}"
	echo "barcode_rev:${barcode_2}"
	echo "length:${length}"
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	if [ -e ./${DIR}/${sample}/${sample}_trimmed_1.fastq.gz ]; then
		echo "already trimmed >>> skip "
	else
		if [ "${direction}" = "single" ]; then
			if [ "${b}" = "barcode" ]; then
				if [ "${length}" = 5 ]; then
					cutadapt -j ${c} -u 29 \
						-o ./${DIR}/${sample}/${sample}_1_tmp.fastq.gz \
						./${DIR}/${sample}/${sample}_1.fastq.gz
					cutadapt -j ${c} -u -102 \
						-o ./${DIR}/${sample}/${sample}_trimmed_1.fastq.gz \
						./${DIR}/${sample}/${sample}_1_tmp.fastq.gz
				elif [ "${length}" = 6 ]; then
					cutadapt -j ${c} -u 30 \
						-o ./${DIR}/${sample}/${sample}_1_tmp.fastq.gz \
						./${DIR}/${sample}/${sample}_1.fastq.gz
					cutadapt -j ${c} -u -101 \
						-o ./${DIR}/${sample}/${sample}_trimmed_1.fastq.gz \
						./${DIR}/${sample}/${sample}_1_tmp.fastq.gz
				elif [ "${length}" = 7 ]; then
					cutadapt -j ${c} -u 31 \
						-o ./${DIR}/${sample}/${sample}_1_tmp.fastq.gz \
						./${DIR}/${sample}/${sample}_1.fastq.gz
					cutadapt -j ${c} -u -100 \
						-o ./${DIR}/${sample}/${sample}_trimmed_1.fastq.gz \
						./${DIR}/${sample}/${sample}_1_tmp.fastq.gz
				elif [ "${length}" = 8 ]; then
					cutadapt -j ${c} -u 32 \
						-o ./${DIR}/${sample}/${sample}_1_tmp.fastq.gz \
						./${DIR}/${sample}/${sample}_1.fastq.gz
					cutadapt -j ${c} -u -99 \
						-o ./${DIR}/${sample}/${sample}_trimmed_1.fastq.gz \
						./${DIR}/${sample}/${sample}_1_tmp.fastq.gz
				fi
			elif [ "${b}" = "default" ]; then
				cutadapt -j ${c} -g "CACCG" \
					-o ./${DIR}/${sample}/${sample}_trimmed_1_tmp.fastq.gz \
					./${DIR}/${sample}/${sample}_1.fastq.gz
				cutadapt -j ${c} -a "GTTTT" \
					-o ./${DIR}/${sample}/${sample}_trimmed_1.fastq.gz \
					./${DIR}/${sample}/${sample}_trimmed_1_tmp.fastq.gz
			fi
			rm ./${DIR}/${sample}/${sample}_*_tmp.fastq.gz
			fastqc -t ${c} --nogroup -o ./${DIR}/${sample}/fastqc/ ./${DIR}/${sample}/${sample}_trimmed_1.fastq.gz
		elif [ "${s}" = "paired" ]; then
			echo "now preparing..."
			#if [ "${b}" = "bar" ]; then
			#	if [ "echo ${#barcode_1}" = 5 ]; then
			#		fwd_1=50
			#		fwd_2=-50
			#		rev_1=50
			#		rev_2=-50
			#	elif [ "echo ${#barcode_1}" = 6 ]; then
			#		fwd_1=50
			#		fwd_2=-50
			#		rev_1=50
			#		rev_2=-50
			#	elif [ "echo ${#barcode_1}" = 7 ]; then
			#		fwd_1=50
			#		fwd_2=-50
			#		rev_1=50
			#		rev_2=-50
			#	elif [ "echo ${#barcode_1}" = 8 ]; then
			#		fwd_1=50
			#		fwd_2=-50
			#		rev_1=50
			#		rev_2=-50
			#	fi
			#	cutadapt -j ${c} -u $((fwd_1)) ./${DIR}/${sample}/${sample}_1.fastq \
			#		> ./${DIR}/${sample}/${sample}_1_tmp.fastq
			#	cutadapt -j ${c} -u $((fwd_2)) ./${DIR}/${sample}/${sample}_1_tmp.fastq \
			#		> ./${DIR}/${sample}/${sample}_trimmed_1.fastq
			#	cutadapt -j ${c} -u ${rev_1} ./${DIR}/${sample}/${sample}_2.fastq \
			#		> ./${DIR}/${sample}/${sample}_2_tmp.fastq
			#	cutadapt -j ${c} -u ${rev_2} ./${DIR}/${sample}/${sample}_2_tmp.fastq \
			#		> ./${DIR}/${sample}/${sample}_trimmed_2.fastq
			#elif [ "${b}" = "def" ]; then
			#	cutadapt -j ${c} -g "CACCG" -a "GTTTT" ./${DIR}/${sample}/${sample}_1.fastq \
			#		> ./${DIR}/${sample}/${sample}_trimmed_1.fastq
			#	cutadapt -j ${c} -g "AAAAC" ./${DIR}/${sample}/${sample}_2.fastq \
			#		> ./${DIR}/${sample}/${sample}_trimmed_2.fastq
			#fi
			#rm ./${DIR}/${sample}/${sample}_*_tmp*.fastq
		fi
	fi
done

echo "~~~~~~~~~~~~trimGuide.sh has been finished~~~~~~~~~~~~~"

exit=0

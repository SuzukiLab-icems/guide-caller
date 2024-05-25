#!/usr/sh
#guide-caller trimGuide program
#Copyright (c) 2023 Yuki NOGUCHI, Jun Suzuki lab
#This code is free software; you can redistribute it and/or modify it
#under the terms of the BSD License (see the file COPYING included with
#the distribution).
#@version: v3.0.0
#@author:  Yuki NOGUCHI
#@contact: jsuzuki AT icems.kyoto-u.ac.jp

DIR=$1		#fastq_directory
sample=$2	#sample
length=$3	#barcode_length
c=$4		#cpu_core

mkdir ${DIR}/${sample}/fastqc

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

fastqc --nogroup -o ./${DIR}/${sample}/fastqc ./${DIR}/${sample}/${sample}_trimmed_1.fastq.gz

exit=0

#!/usr/zsh
echo "guide-caller countGuide program"
echo "Copyright (c) 2023 Yuki NOGUCHI, Jun Suzuki lab"
echo "This code is free software; you can redistribute it and/or modify it"
echo "under the terms of the BSD License (see the file COPYING included with"
echo "the distribution)."
echo "@version: v2.0.0"
echo "@author:  Yuki NOGUCHI"
echo "@contact: jsuzuki AT icems.kyoto-u.ac.jp"

DIR=$1
i=$2
s=$3
VERSION="v2.0.0"
FILE_PATH=~/package/guide-caller/${VERSION}

echo "~~~~~~~~~~~~Welcome to countGuide.sh~~~~~~~~~~~~~"

tail -n+2 ${i} | while read line; do
	sample=$(awk -F ',' '{print $2}' <<< "$line")
	alignment=$(awk -F ',' '{print $3}' <<< "$line")
	alignment="${alignment/#\~/$HOME}"
	direction=$(awk -F ',' '{print $4}' <<< "$line")
	echo "Start analyzing your sample"
	echo "sample:${sample}"
	echo "alignment file:${alignment}"
	echo "length:${direction}"
	echo "FILE_PATH:${FILE_PATH}"
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	if [ -e ${DIR}/${sample}/result/${sample}_.count.txt ]; then
		echo "already mapped >>> skip "
	else
		if [ "${direction}" = "single" ];then
			mageck count \
				-l ${alignment} \
				-n ${sample}_ \
				--fastq ./${DIR}/${sample}/${sample}_trimmed_1.fastq.gz
		elif [ "${direction}" = "paired" ];then
			mageck count \
				-l ${alignment} \
				-n ${sample}_ \
				--fastq ./${DIR}/${sample}/${sample}_trimmed_1.fastq.gz \
				--fastq-2 ./${DIR}/${sample}/${sample}_trimmed_2.fastq.gz \
				--count-pair True
		fi
		mv ${sample}_* ./${DIR}/${sample}/result/
	fi
	python ${FILE_PATH}/src/matrix_shaper.py -o ${DIR}/${sample}/result/${sample} \
		< ${DIR}/${sample}/result/${sample}_.count.txt
done

echo "~~~~~~~~~~~~countGuide.sh has been finished~~~~~~~~~~~~~"

exit=0

#!/usr/sh
#guide-caller countGuide program
#Copyright (c) 2023 Yuki NOGUCHI, Jun Suzuki lab
#This code is free software; you can redistribute it and/or modify it
#under the terms of the BSD License (see the file COPYING included with
#the distribution).
#@version: v2.0.0
#@author:  Yuki NOGUCHI
#@contact: jsuzuki AT icems.kyoto-u.ac.jp

DIR=$1			#fastq_directory
sample=$2		#sample
alignment=$3	#PATH to alignment file
src_path=$4

mkdir ./${DIR}/${sample}/mageck_result/

if [ -e ${DIR}/${sample}/result/${sample}_.count.txt ]; then
	echo "already mapped >>> skip "
else
	mageck count -l ${alignment} -n ${sample} --fastq ./${DIR}/${sample}/${sample}_trimmed_1.fastq.gz
	mv ${sample}* ./${DIR}/${sample}/result/
fi
python ${src_path}/src/matrix_shaper.py -o ${DIR}/${sample}/result/${sample} \
		< ${DIR}/${sample}/result/${sample}.count.txt

exit=0

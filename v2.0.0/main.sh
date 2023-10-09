#!/usr/zsh
echo "guide-caller main program"
echo "Copyright (c) 2023 Yuki NOGUCHI, Jun Suzuki lab"
echo "This code is free software; you can redistribute it and/or modify it"
echo "under the terms of the BSD License (see the file COPYING included with"
echo "the distribution)."
echo "@version: v2.0.0"
echo "@author:  Yuki NOGUCHI"
echo "@contact: jsuzuki AT icems.kyoto-u.ac.jp"
echo "CMD: sh guide-caller/v2.0.0/main.sh -i 231008_crispr_demo.csv -b barcode -e 231008_demo -c 7"
#Opt parse
function usage {
    cat <<EOM
Usage: $(basename "$0") [OPTION]...
    -h          Display help
    -i VALUE	INPUT_info (sample_info.csv)
    -b VALUE	Trimming Start Point (default or barcode): default: Fwd:CACCG, Rev: AAAAC *If you want to use barcode, you put "bar"
    -e VALUE    Experimental_ID (ex. "HN00171248" etc...)
    -c VALUE	Numbers of core
EOM
    exit 2
}

while getopts ":i:b:e:c:s:h" optKey; do
    case "$optKey" in
        i)
          i=${OPTARG}
          ;;
        b)
          b=${OPTARG}
          ;;
		e)
          e=${OPTARG}
          ;;
		c)
          c=${OPTARG}
          ;;
        '-h'|'--help'|* )
          usage
          ;;
    esac
done

#Pre-prep
VERSION="v2.0.0"
FILE_PATH=~/package/guide-caller/${VERSION}
DIR=${e}
sample_ids=$(awk -F "," 'NR>=2 {print $2}' ${i})
direction=$(awk -F ',' 'NR>=2 {print $4}' ${i})
for id in ${sample_ids};do mkdir -p ./${DIR}/${id}/fastqc ./${DIR}/${id}/result;done
#Trimming > QC
sh ${FILE_PATH}/src/trimGuide.sh ${DIR} ${i} ${b} ${c}
#Count
sh ${FILE_PATH}/src/countGuide.sh ${DIR} ${i} ${s}
echo "~~~~ ;-) Bye! ~~~~"
exit 0

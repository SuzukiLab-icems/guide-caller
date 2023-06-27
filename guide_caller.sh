#!/usr/sh

##########################################################################################

##GuideCalling.sh
#Created by NOGUCHI
##########################################################################################
#Defining option
function usage {
    cat <<EOM
Usage: $(basename "$0") [OPTION]...
    -h          Display help
    -i VALUE    ID of experiment (ex. "HN00171248" etc...)
    -f VALUE    Name of Alignment File (ex.Human_GeCKOv2_Library.csv)
    -c VALUE	Numbers of core
EOM
    exit 2
}

#Definition
while getopts ":i:f:c:h" optKey; do
    case "$optKey" in
        i)
          i=${OPTARG}
          ;;
        f)
          f=${OPTARG}
          ;;
		c)
          c=${OPTARG}
          ;;
        '-h'|'--help'|* )
          usage
          ;;
    esac
done

#Environment
LIST=$(ls ./${i})
cp ./alignment_files/${f} ./${i}/
cp matrix.py ./${i}/
cd ./${i}


#FastQC(1st)>Trimming>FastQC(2nd)>Alignment>Count
for SampleID in `echo ${LIST}`
    do
    mkdir ./${SampleID}/${SampleID}_fastqc \
          ./${SampleID}/${SampleID}_fastqc/1st \
          ./${SampleID}/${SampleID}_fastqc/2nd
    fastqc --nogroup -o ./${SampleID}/${SampleID}_fastqc/1st ./${SampleID}/${SampleID}.fastq
    FILE=$(cat ./${SampleID}/${SampleID}.fastq \
        | head -n 2 \
        | tail -n 1 \
        | cut -c 1-5
    )
    if [ "${FILE}" = "ATCGA" ]; then
        cutadapt \
            -j ${c} \
            -u 30 \
            ./${SampleID}/${SampleID}.fastq > 5-trimmed_${SampleID}.fastq
        cutadapt \
            -j ${c} \
            -u -1 \
            5-trimmed_${SampleID}.fastq > trimmed_${SampleID}.fastq
    
    elif [ "${FILE}" = "TAGCT" ]; then
        cutadapt \
            -j ${c} \
            -u 29 \
            ./${SampleID}/${SampleID}.fastq > 5-trimmed_${SampleID}.fastq
        cutadapt \
            -j ${c} \
            -u -2 \
            5-trimmed_${SampleID}.fastq > trimmed_${SampleID}.fastq
    
    elif [ "${FILE}" = "CAAGT" ]; then
        cutadapt \
            -j ${c} \
            -u 28 \
            ./${SampleID}/${SampleID}.fastq > 5-trimmed_${SampleID}.fastq
        cutadapt \
            -j ${c} \
            -u -3 \
            5-trimmed_${SampleID}.fastq > trimmed_${SampleID}.fastq
    
    else
        cutadapt \
            -j ${c} \
            -u 31 \
            ./${SampleID}/${SampleID}.fastq > trimmed_${SampleID}.fastq

    fi
    fastqc --nogroup -o ./${SampleID}/${SampleID}_fastqc/2nd trimmed_${SampleID}.fastq
    mkdir ./${SampleID}/mageck_result
    mageck count \
        -l $f \
        -n ${SampleID} \
        --fastq trimmed_${SampleID}.fastq

    mv ${SampleID}.log ${SampleID}.*.txt *.R *.Rnw ${SampleID}.count_report.Rmd ./${SampleID}/mageck_result
    mv trimmed_${SampleID}.fastq 5-trimmed_${SampleID}.fastq ./${SampleID}/
    
    python matrix.py ${SampleID}/mageck_result ${SampleID}
    

done

exit=0

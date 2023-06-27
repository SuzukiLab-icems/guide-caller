# guide-caller (for the GeCKO v2 sgRNA library analysis pipeline)
This script, guide-caller.sh is an analysis pipeline for Gecko-v2 sgRNA library constructed with FastQC, Cutadapt, and MAGeCK. 
Briefly, the 20bp sgRNA sequences were trimmed from 51bp sequenced read, including plasmid- and adaptor-derived sequences. 
The trimming process was performed twice.
  1st trimming: cutadapt -j 512 -u {N4:28 N5:29 N6:30 N7:31} {Sample}.fastq > {Sample}_(1st)trimmed.fastq
  2nd trimming: cutadapt -j 512 -u {N4:-3 N5:-2 N6:-1 N7:0} {Sample}_(1st)trimmed.fastq > {Sample}_(2nd)trimmed.fastq).
Extracted sgRNAs were mapped to an annotation list by MAGeCK.
In the GeCKO v2 sgRNA library, several sgRNAs potentially annotate multiple genes, mainly family genes (ex., Il11ra family), at the same time.
Therefore, the annotation list was modified to correspond to the sgRNA sequence and all annotated genes.

Reference):
1. FastQC): https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
2. Cutadapt): Martin, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal. 2011; 17: 10-12. 10.14806/ej.17.1.200
3. Mageck): Li, W., Xu, H., Xiao, T., Cong, L., Love, M.I., Zhang, F., Irizarry, R.A., Liu, J.S., Brown, M.A., and Liu, X.S. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology. 2014; 15: 554. 10.1186/s13059-014-0554-4

#!/bin/bash
. env_parallel.bash

# Print a help message
usage() {                                 # Function: Print a help message.
    echo "Usage: $0 [ -b BAM File ] [ -o OUTPUT DIRECTORY ]" 1>&2
    }

# Exit errors
exit_abnormal() {
    usage
    exit 1
}

# Get arguments from the commandline
while getopts b:o: flag
do
    case "${flag}" in
	b)
	    BAM=${OPTARG}
	    ;;
	o)
	    OUTDIR=${OPTARG}
	    ;;
	:)
	    echo "Error: -${OPTARG} requires an argument"
	    exit_abnormal
	    ;;
	*)
	    echo "Error: Unknown  argument"
	    exit_abnormal
    esac
done

if [ $OPTIND -ne 5 ]; then
    echo "Incorrect arguments passed";
    exit_abnormal
fi

NUMCORES=36

echo "Bam: $BAM";
echo "Out Dir: $OUTDIR";

# We need to open many files, so temporarily increase the limit
ulimit -n 4096

# Generate GFF file and spit out cell IDs
mkdir $OUTDIR
mkdir $OUTDIR/individualBams

# Filter bam to remove duplicates to reduce file size
samtools index -@ $NUMCORES $BAM
samtools view -@ $NUMCORES -q 30 -f 3 -bF 3072 $BAM > $BAM.filtered
samtools index -@ $NUMCORES $BAM.filtered
echo "Filtered BAM for quality and removed duplicates"

# Create list of cells to filter based on the metadata and bins and split
Rscript /home/bkw2118/DLP-Overlap-Pipeline/EAC_Sample_Pipeline/CellIDGenerator.R $BAM.filtered $OUTDIR/CellIDList.txt
sinto filterbarcodes -b $BAM.filtered -c  $OUTDIR/CellIDList.txt -p 36 --outdir $OUTDIR/individualBams
echo "Created individual Bam Files"

sort_ind_bam() {
    samtools sort -m 10G -T ${2} -n ${1} | \
	samtools fixmate -rm - - | \
	samtools sort -m 10G -T ${2} - | \
	samtools markdup -r - - | \
	samtools sort -m 10G -T ${2} -n - | \
	bedtools bamtobed -bedpe -mate1 -i - > ${1}.bedpe
}
 
# gnu parallel seems to not recognize that our server is 1 thread per cpu, so set cpu limit to half what you expect 
env_parallel --progress --jobs $(($NUMCORES/2)) sort_ind_bam ::: $(ls $OUTDIR/individualBams/*.bam) ::: $OUTDIR
echo "Filtered Individual BAMs"

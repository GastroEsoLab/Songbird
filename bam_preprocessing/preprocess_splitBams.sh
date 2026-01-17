#!/bin/bash
. env_parallel.bash

# Print a help message
usage() {                                 # Function: Print a help message.
    echo "Usage: $0 [ -b BAM DIRECTORY ] [-t TEMP DIRECTORY]" 1>&2
    }

# Exit errors
exit_abnormal() {
    usage
    exit 1
}

# Get arguments from the commandline
while getopts b:t: flag
do
    case "${flag}" in
	b)
	    BAMDIR=${OPTARG}
	    ;;
	t)
	    TEMPDIR=${OPTARG}
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

NUMCORES=48

echo "Bam Dir: $BAMDIR";

sort_ind_bam() {
    samtools index ${1}
    samtools view -h -q 30 -f 3 -bF 3072 ${1} | \
  samtools sort -m 10G -T ${2} -n - | \
	samtools fixmate -rm - - | \
	samtools sort -m 10G -T ${2} - | \
	samtools markdup -r - - | \
	samtools sort -m 10G -T ${2} -n - | \
	bedtools bamtobed -bedpe -mate1 -i - > ${1}.bedpe

    # Compress for IO when reading in with Songbird
    bgzip -f ${1}.bedpe
}

# gnu parallel seems to not recognize that our server is 1 thread per cpu, so set cpu limit to half what you expect
env_parallel --progress --jobs $(($NUMCORES/2)) sort_ind_bam ::: $(ls $BAMDIR*.bam) ::: $TEMPDIR
echo "Filtered Individual BAMs"

rm ${BAMDIR}*.bam.bai
echo "Removed BAM Index Files"

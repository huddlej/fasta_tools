#!/bin/bash

# Set default contig file to hg19.
sam_extension=".map.gz"
contig_file=~psudmant/genomes/contigs/hg19_contigs.txt

# Get options from the user.
while getopts :c:e: OPTION
do
  case $OPTION in
    c)
      contig_file=$OPTARG
      ;;
    e)
      sam_extension=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac

  shift $OPTIND-1
done

if [[ "$#" -ne "3" ]]
then
    echo "Usage: $0 [-c contigs.txt -s \".map.gz\"] bam_name sam_directory output_dir"
    exit 1
fi

# Get variables.
bam_name=$1
sam_dir=$2
output_dir=$3

# Find all compressed SAM files in the given path, convert them to BAM, and
# sort.
find ${sam_dir} -name "*${sam_extension}" -exec zcat {} \; \
    | samtools view -b -t ${contig_file} - -S | samtools sort - ${output_dir}/${bam_name}.sorted

# Index the resulting BAM file.
samtools index ${output_dir}/${bam_name}.sorted.bam

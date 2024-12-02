#!/bin/bash

# compiled_Data directory  = add the compressed read data here (in .tar)
# reference_genome directory  = has the reference genome in .fa (saaCer3.fa and ty5_6p.fa)
# src directory = will include this automation.sh script

# AT THE END THERE WILL BE TWO NEW DIRECTORIES FOR RESULTS
# COMPILED_IGV_RESULTS = will have the bam file (and .csi index) of each group for use in IGV
# COMPILED_TRANSPOSON_ORIENTATION = will have BED file of all the common singleton reads from each group with BOTH strand information

mkdir -p COMPILED_IGV_RESULTS
mkdir -p COMPILED_TRANSPOSON_ORIENTATION

# STEP1: INDEXING YEAST AND TRANSPOSON GENOME

mkdir -p reference_data/genome_yeast
mkdir -p reference_data/genome_transposon

bowtie2-build reference_data/sacCer3.fa reference_data/genome_yeast/sacCer3
bowtie2-build reference_data/ty5_6p.fa reference_data/genome_transposon/ty5_6p

i=1
# for each group
for file in compiled_data/*.tar
    do 
        echo GROUP ${i}: WORKING ON  $file 

        mkdir Group${i}_file
        mkdir Group${i}_file/data

        mv $file Group${i}_file/data

        cd Group${i}_file
        tar -xvf data/*.tar

        mv *fq.gz data/
        gunzip data/*.gz

        # FASTQC
        mkdir -p results/fastqc
        fastqc data/*.fq -o results/fastqc


        # ALIGN READS TO REFERENCE
        mkdir -p results/bam results/sam
        echo aligning reads to reference for group ${i}

        export BOWTIE2_INDEXES=$(pwd)/../reference_data/genome_yeast
        bowtie2 -x sacCer3 \
            -p 4\
            -1 data/*_1.fq \
            -2 data/*_2.fq \
            -S results/sam/Group${i}_aligned_to_yeast.sam

        export BOWTIE2_INDEXES=$(pwd)/../reference_data/genome_transposon
        bowtie2 -x ty5_6p \
            -p 4\
            -1 data/*_1.fq \
            -2 data/*_2.fq \
            -S results/sam/Group${i}_aligned_to_transposon.sam

        # FILTERING FOR POTENTIAL SINGLETONS
        samtools view -f 9 -F 4 -S -b results/sam/Group${i}_aligned_to_yeast.sam | samtools sort -o results/bam/Group${i}_aligned_to_yeast_singletons.bam
        samtools view -f 9 -F 4 -S -b results/sam/Group${i}_aligned_to_transposon.sam | samtools sort -o results/bam/Group${i}_aligned_to_transposon_singletons.bam

        # GETTING COMMON IDENTIFIERS
        samtools view results/bam/Group${i}_aligned_to_transposon_singletons.bam | cut -f1 | sort -u > results/transposon_singleton_identifiers.txt

        samtools view results/bam/Group${i}_aligned_to_yeast_singletons.bam | cut -f1 | sort -u > results/yeast_singleton_identifiers.txt

        grep -F -x -f results/transposon_singleton_identifiers.txt results/yeast_singleton_identifiers.txt > results/common_identifiers.txt

        rm results/transposon_singleton_identifiers.txt results/yeast_singleton_identifiers.txt

        # GETTING BED FILE OF READS ALIGNED TO TRANSPOSON / YEAST
        mkdir -p results/bed

        bedtools bamtobed -i results/bam/Group${i}_aligned_to_yeast_singletons.bam > results/bed/Group${i}_aligned_to_yeast_singletons.bed
        bedtools bamtobed -i results/bam/Group${i}_aligned_to_transposon_singletons.bam > results/bed/Group${i}_aligned_to_transposon_singletons.bed

        grep results/bed/Group${i}_aligned_to_transposon_singletons.bed -F -f results/common_identifiers.txt > results/bed/Group${i}_aligned_to_transposon_common_identifiers.bed
        grep results/bed/Group${i}_aligned_to_yeast_singletons.bed -F -f results/common_identifiers.txt > results/bed/Group${i}_aligned_to_yeast_common_identifiers.bed

        # OBTAINING BED FILE CONTAINING GENOME COVERAGE INFORMATION
        samtools faidx ../reference_data/sacCer3.fa

        cut -f1,2,3 results/bed/Group${i}_aligned_to_yeast_common_identifiers.bed > temp.bed

        bedtools genomecov -i temp.bed -g ../reference_data/sacCer3.fa.fai -bg > results/bed/genome_coverage.bed

        rm temp.bed

        # CREATION OF BAM FILE FOR VIEWING ON IGV
        mkdir -p results/for_IGV
        samtools view -h -b -L results/bed/genome_coverage.bed results/bam/Group${i}_aligned_to_yeast_singletons.bam | samtools sort -o results/for_IGV/GROUP${i}_subset_of_reads_aligned_to_yeast_forIGV.bam --write-index

        # INFERRING THE ORIENTATION OF THE TRANSPOSON INSERTION AT EACH INSERTION SITE
        sort results/bed/Group${i}_aligned_to_yeast_common_identifiers.bed -k4,4 > results/temp_yeast_reads_info_identifier_sorted.txt
        # Sorting the reads aligned to transposon according to identifier names to match order of reads in the file above, and extracting the "strand" field into a temporary file
        sort results/bed/Group${i}_aligned_to_transposon_common_identifiers.bed -k4,4 | cut -f6 > results/temp_transposon_reads_info_identifier_sorted.txt

        # Combine the file containing information of each read aligned to yeast genome and the
        # additional file containing just the transposon strand which the read with the identical identifier mapped to, for each row
        # Also, sorting of reads according to chromosome, then base pair position
        paste results/temp_yeast_reads_info_identifier_sorted.txt results/temp_transposon_reads_info_identifier_sorted.txt | \
            sort -k1,1 -k2,2n > results/GROUP${i}_combined_reads_information_for_orientation_determination.txt

        mv results/GROUP${i}_combined_reads_information_for_orientation_determination.txt ../COMPILED_TRANSPOSON_ORIENTATION

        rm results/temp_yeast_reads_info_identifier_sorted.txt results/temp_transposon_reads_info_identifier_sorted.txt


        # Move results to final output directory
        mv results/for_IGV/GROUP${i}* ../COMPILED_IGV_RESULTS
        cd ..
        echo COMPLETED FOR GROUP ${i}
        let i=i+1
    done







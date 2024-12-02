Developed and implemented a fully automated Bash pipeline for analyzing high-throughput sequencing data to investigate transposon insertion orientation and genomic locations in yeast.

Automation.sh workflow:

Index yeast (sacCer3) and transposon reference genomes using bowtie2-build.
Align paired-end sequencing reads to yeast and transposon references with bowtie2.
Access read quality of raw FASTQ files via FastQC.
Use samtools to filter singleton reads aligning uniquely to yeast and transposon genomes for transposon-specific analysis.

Convert aligned reads to BED format using bedtools and generated genome coverage files for visualization and interpretation.
Generate IGV-compatible BAM files with indexed genome coverage for seamless inspection of alignment results.
Combine yeast and transposon strand information to determine the orientation of transposon insertions at genomic integration sites.

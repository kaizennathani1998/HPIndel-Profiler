# HPIndel-Profiler
Command-line tool for detecting and quantifying indel events in homopolymer and non-homopolymer regions from BAM files using BED-defined genomic intervals. Outputs a .txt file with tab delimated information of input file, hp indel rates, non-hp indel rates and ratio

# Inputs required:
1) Reference genome fasta and fasta index file (Eg: hg19.fa, hg19.fai.fa)
2) Homopolymer region BED file
3) Homopolymer region BED file

# Usage
./hp_indelrate.sh <BAM_FILE> <resources/HP_region.bed> <resources/NON-HP_region.bed> <output_file_path>

### Note: Input BED files should be in resources folder only 

# Apptainer run
apptainer pull hp_indelrate.sif docker://docker.io/kaizennathani1998/hp-indelrate:latest

apptainer run \\
  hp_indelrate.sif \\
  <BAM_FILE> \\
  <HP_BED> <NON-HP_BED> \\
  <output summary fille path.txt>



<img width="1144" height="763" alt="image" src="https://github.com/user-attachments/assets/0a3af9e0-1dbb-4547-abc5-bd286d3bacfc" />

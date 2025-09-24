#!/bin/bash

# HackBio Internship: Stage 0 
# Chloe Nichole Calica
# August 30, 2025

# Project 1: BASh Basic

# 1. Print your name
echo ‘Chloe Nichole Calica’

# 2. Create a folder titled your name
mkdir Chloe

# 3. Create another new directory titled biocomputing and change to that directory with one line of command
mkdir biocomputing && cd biocomputing

# 4. Download 3 files
curl -O https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna -O https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk -o wildtype2.gbk https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk

# 5. Move the .fna file to the folder titled your name
mv wildtype.fna ../Chloe

# 6. Delete the duplicate gbk file
rm wildtype2.gbk

# The .fna file is actually from a bacteria, and it should definitely have a TATA (tata) box for initiating gene transcription. A molecular biologist is trying to understand the implication of dual TATA sequences. The files got mixed up and we are not sure which is wildtype and which is mutant. The mutant should have “tatatata” while the normal should have just “tata”. Can you confirm if the file is mutant or wild type

# 7. Confirm if the .fna file is mutant or wild type (tatatata vs tata)
if grep -q ‘tatatata’ wildtype.fna; then
	echo ‘mutant’
elif grep -q ‘tata’ wildtype.fna; then
	echo ‘wildtype’
else
	echo ‘No TATA box or dual TATA present.’
fi

# 8. If mutant, print all matching lines into a new file
grep 'tatatata' ../Chloe/wildtype.fna > dual_tata.txt

# 9. Count number of lines (excluding header) in the .gbk file

# Header includes everything before FEATURES so to count excluding the header, we find the line that contains the first instance of FEATURES then we count all the lines starting from there

awk '/FEATURES/ {found=1} found' wildtype.gbk | wc -l

# 10. Print the sequence length of the .gbk file. (Use the LOCUS tag in the first line)

# the sequence length is the third field in the LOCUS line
awk '/LOCUS/ {print $3}' wildtype.gbk

# 11. Print the source organism of the .gbk file. (Use the SOURCE tag in the first line)

# sed is the stream editor and works line by line
# arg -n tells sed to not print every line automatically, just the matches
# s/…/…/, substitute command: s/pattern/replacement
# ^SOURCE, matches lines that start with SOURCE
# \s, whitespace characters
# \+, one or more
# //, replacement is empty so removes SOURCE and the spaces after it
# After substitution, we are left with just the organism name

sed -n 's/^SOURCE\s\+//p' wildtype.gbk

# 12. List all the gene names of the .gbk file. Hint {grep '/gene='}

# use grep to find all the lines with genes
# use awk to clean the lines and only keep gene names
# arguments -F to specify separators, in this case the quotation marks which leaves the gene name in the 2nd field

grep '/gene=' wildtype.gbk | awk -F'"' '{print $2}'

# 13. Clear your terminal space and print all commands used today
clear && history

14. List the files in the two folders
# ~, home directory
# -R, recursively list all subdirectories

ls -R ~

# Project 2: Installing Bioinformatics Software on the Terminal

# 1. Activate your base conda environment

# Download the Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer
bash Miniconda3-latest-Linux-x86_64.sh

# Add Conda to my PATH for this current session
export PATH="$HOME/miniconda3/bin:$PATH"

conda init
conda version

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge


# 2. Create a conda environment named funtools
conda create -n funtools python=3.10

# 3. Activate the funtools environment
conda activate funtools

# 4. Install Figlet using conda or apt-get
sudo apt install figlet

# 5. Run figlet <your name>
figlet chloe

# 6. Install bwa through the bioconda channel
conda install -c bioconda bwa

# 7. Install blast through the bioconda channel
conda install -c bioconda blast

# 8. Install samtools through the bioconda channel
conda install -c bioconda samtools

# 9. Install bedtools through the bioconda channel
conda install -c bioconda bedtools

# 10. Install spades.py through the bioconda channel
conda install -c bioconda spades.py

#11. Install bcftools through the bioconda channel
conda install -c bioconda bcftools

#12. Install fastp through the bioconda channel
conda install -c bioconda fastp

#13. Install multiqc through the bioconda channel
conda install -c bioconda multiqc

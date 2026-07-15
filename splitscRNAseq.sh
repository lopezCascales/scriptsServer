curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
info: downloading installer

Welcome to Rust!

This will download and install the official compiler for the Rust
programming language, and its package manager, Cargo.

Rustup metadata and toolchains will be installed into the Rustup
home directory, located at:

  /home/lopez_cascales/.rustup

This can be modified with the RUSTUP_HOME environment variable.

The Cargo home directory is located at:

  /home/lopez_cascales/.cargo

This can be modified with the CARGO_HOME environment variable.

The cargo, rustc, rustup and other commands will be added to
Cargo's bin directory, located at:

  /home/lopez_cascales/.cargo/bin

This path will then be added to your PATH environment variable by
modifying the profile files located at:

  /home/lopez_cascales/.profile
  /home/lopez_cascales/.bashrc

You can uninstall at any time with rustup self uninstall and
these changes will be reverted.

Current installation options:


   default host triple: x86_64-unknown-linux-gnu
     default toolchain: stable (default)
               profile: default
  modify PATH variable: yes

1) Proceed with installation (default)
2) Customize installation
3) Cancel installation


- rustc 1.76.0 (07dca489a 2024-02-04)


Rust is installed now. Great!

To get started you may need to restart your current shell.
This would reload your PATH environment variable to include
Cargo's bin directory ($HOME/.cargo/bin).

To configure your current shell, run:
source "$HOME/.cargo/env"




####################################################################################
#! /bin/env bash

## Description: Split scRNAseq Fastq by barcodes 
## Date: May 14, 2023
## Version: 1.0
## Contributors: Mayte Lopez-Cascales
## Input: FASTQ files and barcodes.txt
## Output: Bam files splits by cell barcodes

# BEFORE RUNNING THIS SCRIPT

## 1. Run -> chmod u+x splitscRNAseq.sh
## 2. Execute script in your folder with -> bash SCRIPTsplitfastqbybarcodes.sh  (./SCRIPTsplitfastqbybarcodes.sh)

experimentID="yourexperiment"
echo "---------ID = $experimentID ---------"

## Output directory
rm -rf output
mkdir -p output/{fastq,bamfiles}


# FASTQ file statistics
seqkit stats *.fastq.gz &> output/${experimentID}_stats.txt

# Run FASTQC
echo "---------FASTQC---------"
fastqc *.fastq.gz -o output/fastqc -q

### merge fastq IN A BIG FASTQ (If this step its neccesary, depends of your experiment)

#cat *_L001_R1.fastq *_L002_R1.fastq > bigfile_R1.fastq
#cat *_L001_R2.fastq *_L002_R2.fastq > bigfile_R2.fastq

# example
# R1
# cat ERA-MC-s159_HW5N3BGXN_S1_L001_R1_001.fastq.gz ERA-MC-s159_HW5N3BGXN_S1_L002_R1_001.fastq.gz ERA-MC-s159_HW5N3BGXN_S1_L003_R1_001.fastq.gz ERA-MC-s159_HW5N3BGXN_S1_L004_R1_001.fastq.gz > s159_R1.fastq.gz

# R2
#cat ERA-MC-s159_HW5N3BGXN_S1_L001_R2_001.fastq.gz ERA-MC-s159_HW5N3BGXN_S1_L002_R2_001.fastq.gz ERA-MC-s159_HW5N3BGXN_S1_L003_R2_001.fastq.gz ERA-MC-s159_HW5N3BGXN_S1_L004_R2_001.fastq.gz > s159_R2.fastq.gz

# Unzip the .gz files
echo "---------FASTQC---------"
gunzip *.fastq.gz


### SPLIT the fastq files by barcodes 

#for i in $(cat barcodes.txt);do seqkit grep -r -p $i s159_R1.fastq > $i_R1.fastq ;done
#for i in $(cat barcodes.txt);do seqkit grep -r -p $i s159_R2.fastq > $i_R2.fastq ;done

for experimentID in $(cat barcodes.txt);do seqkit grep -r -p ${experimentID} *_R1.fastq > ${experimentID}_R1.fastq ;done
for experimentID in $(cat barcodes.txt);do seqkit grep -r -p ${experimentID} *_R2.fastq > ${experimentID}_R2.fastq ;done

### STAR alignment (Generate Bam files)




### 

#################################################################

### How to split bam files by strand (only if you need)
# https://www.biostars.org/p/92935/
# https://www.biostars.org/p/196746/#196755
# https://www.biostars.org/p/179035/

#for file in ./*.bam
#do
#echo $file 
#samtools view -b -F 16 $file > ${file/.bam/_rev.bam}
#done

#for file in ./*.bam
#do
#echo $file 
#samtools index $file
#done

#for file in ./*.bam
#do
#echo $file 
#samtools view -b -f 16 $file > ${file/.bam/_fwd.bam}
#done

#for file in ./*fwd.bam
#do
#echo $file 
#samtools index $file
#done

#####################################

# filter by quality

for file in ./.*.bam;
do echo $file
samtools view -b -q 10 $file > ${file/.bam/filtered.bam}; 
done

#####################################

# sort by position

for file in ./.*filtered.bam; 
do 
samtools sort $file > ${file/.filtered.bam/sortedPos.bam};
done

#####################################

# samtools index

for file in ./.*SortedPos.bam; 
do 
samtools index $file;
done

#####################################


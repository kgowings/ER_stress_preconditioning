#The High Throughput Genomics Core sequenced my samples at the Huntsman Cancer Institute at the Univesity of Utah (https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg) and I was provided with compressed raw FASTQ files 
#I performed the following code in the University of Utah's HCP clusters (https://www.chpc.utah.edu/resources/HPC_Clusters.php) 

#SOFTWARE- All software needed for this project was already downloaded on the clusters, and I just needed to load the required modules  
    #to make sure that it is downloaded run:
      module spider SOFTWARE
    #this will provide you with information on whether it is downloaded and how to run it

#SLURM JOBS: You can't run any computationally heavy jobs on your home node.
#Your home node is primarily just used for making folders, moving files, deleting files, and whatnot. For computationally intensive jobs, you will need to submit to sbatch using slurm - can't run it on the home node - need to submit to the compute node.
#A lot more information about slurm and how to use it on chpc website https://www.chpc.utah.edu/documentation/software/slurm.php
#EXAMPLE of a SLURM script 

#to create a job, you will need to create a script for SLURM using vim
vim JOB.sh 
#once in vim, press "i" to enter insert mode and type in (must include the # signs in front of the SLURM description) 
#To view a list of accounts and partitions that are available to you, run command: myallocation

#!/bin/tcsh
#SBATCH --partition=PARTITION_NAME		
#SBATCH --account=ACCOUNT_NAME					
#SBATCH --nodes=2							              	#this number is variable - I am cautious and choose a small number - don't wanna get kicked off 
#SBATCH -o %j-%N.out					                #this will provide you a file (job#.out) with the normal output from the job
#SBATCH -e %j-%N.err					              	#this will provide you a file (job#.err) with any errors produced from the job
#SBATCH -J Job Name   				          	    #name your job - makes finding it in the squeue easier
#SBATCH --time=6:00:00						      	    #you don't need but can choose to include a time limit 
#SBATCH --mail-type=FAIL,ENDjkh,ljhl			    #telling the server to email you if the job fails
#SBATCH --mail-user=EMAIL@utah.com		#telling the server where to email you
#	
#*All of the pound (#) signs before !/bin/tcsh and before SBATCH is all needed in your script - this tells SLURM HOW to run the job even though it is not PART of the job*

#This is where you would type out your job- I will only be showing this section of the script in the following steps 
#When you submit a job, the environment it is in is completely different than your home node, so make sure that you load any modules that you will use in your job script here as well as defining the path of the files that you are performing the job as well as the path where you want the destination files

#*Note: you can perform multiple steps in one script. Ex: you can trim, then index, then align all in one script so you only have to submit one job. However, you might want to start with a script per job, step-by-step, in case you get any errors you know exactly what it was from***

#to exit vim, hit "esc" then type ":wq" to save and exit

#SUBMITTING A JOB: Once you have created your job script you need to run it with sbatch
sbatch run.sh

#CHECKING ON YOUR JOB: you can check on your job after you have submitted it, to see if it is running or pending, and how long it has been running
squeue -u UNID

#SCRATCH FOLDER: make scratch file for fastq files (and all future files - there is a lot more storage there - but unsused files will get deleted after 60 days -so make sure you have a backup somewhere else!)
mkdir DIRECTORY_LOCATION/FASTQ
cd DIRECTORY_LOCATION/FASTQ

#DOWNLOADING FASTQ Files from GNOMEX Website using FDT Command Line
#if files are larger than 1GB (which they will be most certainly) it can be faster/easier to use FDT Command Line
#First, download the fdtCommandLine.jar app from http://hci-bio-app.hci.utah.edu/fdt/
#then put it in your home directory: ~/fdtCommanLine.jar
#the get the location of the files you want through the genomex website
#execute this command with your own path variables:

java -jar ~/fdtCommandLine.jar -noupdates -pull -r -c hci-bio-app.hci.utah.edu -d DIRECTORY_LOCATION/FASTQ GENOMEX_LOCATION/folder  

#Example
#java -jar ~/fdtCommandLine.jar -noupdates -pull -r -c hci-bio-app.hci.utah.edu -d /scratch/kingspeak/serial/UNID/FASTQ /scratch/fdtswap/fdt_sandbox_gnomex/82b224a4-e543-4161-8e0b-2a9e40f43222/PROJECT_ID

###STEP 1: EXTRACT gzipped files
#need seqtk software 
#this step is quick and easy enough that you don't need to submit a job through SLURM, can do it on your home node

gzip -d -c filename.fastq.gz > filename.fastq

	#Example
	#gzip -d -c RAL737_noHS_L001_R1_001.fastq.gz > RAL737_noHS_L001_R1.fastq
	#gzip -d -c RAL737_noHS_L002_R2_001.fastq.gz > RAL737_noHS_L002_R2.fastq

###STEP 2:CONCATENATE lanes 
#In this case, each sample was sequenced in two separate lanes, L001 and L002, of the same flow cell. 
#There are two ways to combine the lanes: concatenate the Fastq files prior to alignment, or align individually, then concatenate the alignments (using Linux samtools merge). There are technical reasons for doing one or the other in certain situations, but in this case we I did the first one for simplicity.

cat filename*.fastq > filename.fastq

	#Example
    #Files being combined: RAL387_HS_L001_R2.fastq & RAL387_HS_L002_R2.fastq
	#cat RAL737_noHS_L*_R1.fastq > RAL737_noHS_R1.fastq


###STEP 3:TRIMMING READS

module load seqtk/040218
#trim fastq files with seqtk 
seqtk trimfq filename.fa > filename.trimmed.fa

    #Example- included SLURM code above
    
    #module load seqtk/040218
    #seqtk trimfq RAL819_noHS_R1.fastq > RAL819_noHS_R1.trimmed.fastq
    #seqtk trimfq RAL819_noHS_R2.fastq > RAL819_noHS_R2.trimmed.fastq

###STEP 4: INDEX the genome

#First you need to upload the genome to your directory
#THEORETICALLY you would only need the transcriptome for RNA seq. However, issues arise when you do that! So make sure that even for RNAseq you get the GENOME not the TRANSCRIPTOME. (Make sure you have the correct organism, as well as the correct build that you want) 
#I used Drosophila melanogaster reference genome (assembly BDGP6.28, Ensembl release 109) - https://ftp.ensembl.org/pub/release-109/fasta/drosophila_melanogaster/dna_index/
#Software needed- bowtie2/2-2.2.9 & python/3.10.3 (some bowtie functions depend on python) 
#For more about Bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

module load bowtie2/2-2.2.9
module load python/3.10.3
bowtie2-build -f DIRECTORY_LOCATION/BDGP6.28_release109/D_melanogaster.BDGP6.32.dna.toplevel.fa DIRECTORY_LOCATION/BDGP6.28_release109/Drosophila_melanogaster.BDGP6.32.dna.toplevel.genome

#your job output should be 6 different files - "bowtie2-build outputs a set of 6 files with suffixes .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2."

###STEP 5: ALIGNMENT 
#This step will take a few hours, even on the server and even if you submit several sbatch jobs 

module load bowtie2/2-2.2.9
module load python/3.10.3
bowtie2 -x indexed genome -1 trimmed fatsq -2 paried trimmed fastq -S output_align.sam

    #EXAMPLE
    #module load bowtie2/2-2.2.9
    #module load python/3.10.3
    #bowtie2 -x /scratch/general/vast/u6004332/BDGP6.28_release109/Drosophila_melanogaster.BDGP6.32.dna.toplevel.genome -1      /scratch/general/vast/u6004332/RawFASTQ_HSDGRP/RAL69_HS_R1.trimmed.fastq -2 /scratch/general/vast/u6004332/RawFASTQ_HSDGRP/RAL69_HS_R2.trimmed.fastq -S /scratch/general/vast/u6004332/SAM/RAL69_HS.sam


###STEP 5: SAMTOOLS 
#convert and sort the sam files into bam files 
#Software- samtools/1.5

module load samtools/1.5
samtools view -Sb output_align.sam > output.bam

    #EXAMPLE
    #module load samtools/1.5
    #samtools view -Sb /scratch/general/vast/u6004332/SAM/RAL195_HS.sam > /scratch/general/vast/u6004332/BAM/RAL195_HS.bam

#Need to sort your bam file that you just created - do this with samtools as well

module load samtools/1.5
samtools sort output.bam -o output.sorted.bam

    #EXAMPLE
    #module load samtools/1.5
    #samtools sort /scratch/general/vast/u6004332/BAM/RAL195_HS.bam -o /scratch/general/vast/u6004332/BAM/RAL195_HS.sorted.bam

#You can also wait to combine lanes and merge the sorted bam files with samtools 
samtools merge output.bam input.bam input.bam

#If jobs keep getting killed due to preemption, try adding
#SBATCH --requeue
    #this will automatically requeue to job to be run if it keeps getting kicked off by other jobs with higher priority
#Can check for jobs with squeue -u UNID but  can also do scontrol show job JOBID for more detailed description of your job (like predicted finish time)
#Cancel job with scancel <JOBID>

***Now that you have your sorted BAM files, you are ready to move onto DEseq2!***



    


















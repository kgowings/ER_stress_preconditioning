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





















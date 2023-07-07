###STEP 1: GWA
#Gather the following 5 files into one DIRECTORY on your server- 
  #Plink formatted genotype files (BED/BIM/FAM) 
    #Download from http://dgrp2.gnets.ncsu.edu/data.html
    #Don't delete anything from the dgrp2.fam file. Enter the phenotype result data from the screen into the last column. Put -9 for lines with no data. 
  #input.txt
    # two columns: strain result
  #GEMMA.txt 


#In server or file manager of your choice
#Open GEMMA.txt in DIRECTORY
  #Change name of input file to match file name 
  #Make sure account and partition are correct
  #Change output file name to name desired
  #Save
#Change GEMMA.txt extension to .sh (can right click and rename in winscp)

#Ready to run!
#On server, in DIRECTORY run 
sbatch GEMMA.sh

#To check status type
Squeue -u unid

#Gather the following 5 files into one DIRECTORY on your server- 
  #Plink formatted genotype files (BED/BIM/FAM) 
    #Download from http://dgrp2.gnets.ncsu.edu/data.html
    #Don't delete anything from the dgrp2.fam file. Enter the phenotype result data from the screen into the last column. Put -9 for lines with no data. 
  #input.txt
    # two columns: strain result
  #GEMMA.txt - .txt sbatch file containing the following: 

/usr/bin/terminal
srun --time=3:00:00 --nodes=2 --account=XXX --partition=XXX --pty /bin/bash -l
cd /DIRECTORY_location
#your path must already have GEMMA (https://github.com/genetics-statistics/GEMMA) 
module load gemma
gemma -bfile dgrp2 -gk 1 -o dgrp2_matrix
mv /DIRECTORY_location/output/dgrp2_matrix_cXX.txt /DIRECTORY_location
cd /DIRECTORY_location
gemma -bfile dgrp2 -k dgrp2_matrix.cXX.txt -maf 0.001 -miss 0.2 -lmm 4 -n 1 -o OUTPUT_FILE_NAME
cd /DIRECTORY_location/output
awk '$12 < 10^-4' OUTPUT_FILE_NAME.assoc.txt > OUTPUT_FILE_NAME_assoc.top.txt

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

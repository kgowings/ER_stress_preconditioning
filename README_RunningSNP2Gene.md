#These can be run on your desktop or the server. Very low computational power is needed. 
#This example was run on a PC desktop using Gitbash 

#Gather the following files in one folder- 
#dgrp.fb557.txt - download at http://dgrp2.gnets.ncsu.edu/data/website/dgrp.fb557.annot.txt
#GWASresults.txt -  OUTPUT_FILE_NAME.assoc.txt output file from https://github.com/kgowings/ER_stress_preconditioning/blob/main/DGRP%20GWA%20(Table%201)
#SNP2gene.py - https://github.com/kgowings/ER_stress_preconditioning/commit/bbcb14cb7431aebbc7316b48135d30b1d26f33c3
#SNP2Gene4GSEA.py - https://github.com/kgowings/ER_stress_preconditioning/blob/main/SNP2Gene4GSEA.py

#Open SNP2gene_2.py in Notepad++ 
line 14: change GWAfile to match file name and location to GWASresults.txt file name and location
line 15: change SNPfile destination to match your folder 
line 16: SNP2gene_OUTPUT, change destination to match your folder and change name to desired output name 
Save 

Open Gitbash to the folder and run the following code: 
#If Python isn't your path, run:
#PATH=$PATH:/c/Python27/
python SNP2gene_2.py

Open SNP2Gene4GSEA.py in Notepad++
line 36-40, edit to match your folder locations 
36: match to GWASresults.txt file name
38: match SNP2gene_OUTPUT.txt file name from previous step 
41&42: name your output files 

Open Gitbash to folder and run following code: 
#If python isn't your path run: PATH=$PATH:/c/Python27/
python SNP2Gene4GSEA_2.py


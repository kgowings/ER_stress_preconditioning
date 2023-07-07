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

###STEP 2: Assign SNPS to fly genes and prepare file for GSEA
#This can be run on your desktop or the server. Very low computational power is needed. 
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

###STEP 3: Run GSEA
###Code originally published: https://gsajournals.figshare.com/articles/dataset/Supplemental_Material_for_Chow_et_al_2019/9808379?file=17599163
###Code created in the Goodman lab and published in association with the following manuscript: https://doi.org/10.1534/g3.119.400722
#These can be run on your desktop or on the server. Very low computational power needed. 
#On desktop, use Gitbash (or Mac equivalent) 

#Need the following files in one directory- 
  #all_droso_genes.txt <- [expression_file]
  #GSEAsourceCode.txt
  #term2id.bp.mf.cc.txt
  #SNP2Gene4GSEA_output_terse.txt <- [annotation_file]

#Fifteen individual Java class files are listed in the associated file “GSEAsourceCode.txt”
#Compile the class files into GSEA.jar using any compiler of your choosing.

#To run GSEA.jar, use the following command:

java -jar GSEA.jar [expression_file] [annotation_file] [output_file]

#from the directory holding the jar file, expression file, and annotation file.  
#The expression_file is the one created from the GWAS output file, consisting of two columns (FBgn# and -log(GWAS variant ID p-value)).  
#The annotation_file is a tab-delimited file that links the Entry_ID (i.e., FBgn#) with the functional groups (i.e., GO ID) to test against.  
#The output_file is where you want the results to go.  I usually call this OUTPUT.txt and rename it later.

#The output format will be a tab-delimited file as follows
#[experiment_name] [functional_group] [p_value] [es_score] [#genes] [genes]

#I open the OUTPUT.txt file in Excel and add these headings in the 1st row.  
#The es_score is the marker of whether the gene list was concentrated at the top (positive) or bottom (negative).  
  #When the input data consists of only one experiment, as is the case with GWAS data, anything with a negative es_score means the gene list was concentrated at the bottom of the list, which corresponds to things not    significantly enriched.  
#The p-value is the significance of the es_score, computed using the dynamic programming method of Keller et al (PMID 17683603).  
#The #genes tells you the number of genes contributing to the maximum es score found and the genes tells you the ids of the genes that are in that list as well (seperated by ";").

Open OUTPUT.txt in excel 
#With the Excel file, you can use the term2id.bp.mf.cc.txt file with VLOOKUP to assign GO category names to the ID#’s. 
#Then I do a series of sortings to delete rows with negative enrichment scores and #genes less than a cutoff value (usually 3-5) but depends on your output data. 
#Finally, I use the last column (genes) and assign gene names to each and the (-log(GWAS variant ID p-value)).  
#I only do this if I want to show the genes found within each enriched GO term.

###STEP 4: GSEA R plot

#Create GSEAplot.txt: 
  #Example
    #Term	pvalue	ES	No_genes
    #male meiosis I	0.015876766	0.518760295	8
    #negative regulation of growth	0.01509331	0.521337929	6 
  #^Sort file however you want terms to be sorted on plot (by p_value or es score) 

Run GSEAplot.R

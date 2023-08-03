###General DESeq2 workflow for differential expression analysis of RNAseq data### 
	
#Helpful website for writing your own DESEQ work flow and undertanding the steps, I highly recommend reading through this even if you use the following code: http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
#If you have FASTQ files and want to align with bowtie2 start by with https://github.com/kgowings/ER_stress_preconditioning/blob/main/RNAseq_DESeq/README_Alignment.md

#Before beginning make sure you have your bam files, reference genome, and sampletables all in one folder in the server ("/DIRECTORY/")
#SampleTable is a csv file containing information about each sample file 
	#First column is the sample name
	#Second column the file (.bam)
	#Remaining columns are sample metadata which will be stored 
	#Example sampleTable.csv
		##sample	fileName	group	treatment 
		##RAL69_HS	RAL69_HS.sorted.bam	beneficial	heatshock
		##RAL819_noHS	RAL819_noHS.sorted.bam	detrimental	control
		## ...
		 
#sampleTable_DESeq.csv is a csv file containing information about the samples being analyzed/compared- needed for step 2 when making the colData(sampleTable)
	#Example of sampleTable_DESeq.csv
	#RAL304_noHS.sorted.bam,,detrimental
	#RAL335_noHS.sorted.bam,detrimental
	#RAL737_noHS.sorted.bam,detrimental
	#RAL819_noHS.sorted.bam,detrimental
	#RAL195_noHS.sorted.bam,detrimental
	#RAL69_noHS.sorted.bam,beneficial
	#RAL93_noHS.sorted.bam,beneficial
	#RAL359_noHS.sorted.bam,beneficial
	#RAL387_noHS.sorted.bam,beneficial
	#RAL409_noHS.sorted.bam,beneficial


#I used FastX3 - https://www.chpc.utah.edu/documentation/software/fastx.php
#Allows you to interact with remote linux systems graphically, which is useful when running RStudio in a server
	
#Launch an xterm session, launch interactive node, then type:
module load R
module load RStudio 
rstudio 

#downloading packages needed
#its best to try to download everything you might need at the beginning so there is no dependency problems 
#you only have to download the programs once, then they should be on your server for future uses
	##To install core packages, type the following in an R command window
	if (!require("BiocManager", quietly = TRUE))
		install.packages("BiocManager")
	BiocManager::install()
	##Install specific packages with
	BiocManager::install(c("GenomicFeatures","Rsamtools", "pheatmap", "RColorBrewer", "PoiClaClu","org.Hs.eg.db", "ReportingTools", "Gviz", "sva", "magrittr"))


###STEP1-LOAD/PREPARE FILES

indir <- file.path("/DIRECTORY/")
#this should be the path to access your sample table, bam files, and gtf reference genome 

list.files(indir)
#Load sample tables
sampleTable_path <- file.path(indir, "sampleTable.csv")
#OR you can use:
sampleTable <- read.csv(sampleTable_path)

sampleTable_DESeq_path <- file.path(indir, "sampleTable_DESeq.csv")

#Define variable- filename
#The following code will change based on the setup of your sample table
filename <- file.path(indir, paste0(sampleTable$fileName, ".bam"))
file.exists(filename)
#it should say TRUE for each bam file you expect to have. If it doesn't you need to re-examine your code before moving forward

#Define variable- gtffile
gtffile <- file.path(indir,"GTF.gtf")
file.exists(gtffile)
#Loads your refernce genome, make sure this is the same reference used in your alignment step
#Since I used Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa for alignment, I used Drosophila_melanogaster.BDGP6.28.102.gtf
#Downloaded gtf from Ensembl: http://useast.ensembl.org/Drosophila_melanogaster/Info/Index
#More information on GTF format: https://ftp.ensembl.org/pub/release-109/gtf/drosophila_melanogaster/README

#To unzip GTF file use
gzip -d -c FILENAME.gtf.gz > FILENAME.gtf


###Step 2- CREATE A COUNT MATRIX
#the next step is computationally heavy so an interactive batch job must be run 
#exit out of RStudio with ctrl+C in xterm 
#start interactive session on the compute node:
salloc --time=2:00:00 --ntasks=2 --nodes=1 --account=XXX --partition=XXX 
module load R
module load RStudio 
rstudio

##Step 2-OPTION 1-SummarizedExperiment Input <- Have used before but I used OPTION 2 for final results in the manuscript
library("Rsamtools")
library("GenomicFeatures")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
ebg <- exonsBy(txdb, by="gene")
library("GenomicAlignments")
#read counting step
se <- summarizeOverlaps(features=ebg, reads=filename, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
	#Counting mode "Union" = those reads which overlap any portion of exactly one feature are counted. 
	#singleEnd =FALSE = experiment produced paired-end reads
	# ignore.strand = TRUE = Protocol was not strand-specific
	#fragments = TRUE = count reads with unmapped pairs. This is only for paired-end experiments.
	#this line of code will change based on your sequencing protocol
	#this step takes a few minutes to run- one of the more computationally heavy steps 
colData(se) <- DataFrame(sampleTable_DESeq)
	#colData slot, so far empty, should contain all the metadata. 
	#This line assigns the sampleTable as the colData of the summarized experiment by converting it into a DataFrame and using the assignment function

#SummarizedExperiment input= http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#summarizedexperiment-input
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ treatment)
dds
#another example of the design:
# ~ cell + sample.name meaning that we want to test for the effect of treatment+infection (included in sample name) controlling for the effect of different cell line 


##Step 2-OPTION 2-Count matrix input<- option I used/*recommend*##
#featureCounts: a ultrafast and accurate read summarization program

library("Rsubread")
	#documentation- https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf
	#input- SAM/BAM files and an annotation file including chromosomal coordinates of features (gtf/gtf.gz)
	#output- numbers of reads assigned to features (or meta-features) & stat info for the overall summrization results
	#feature- each entry in the provided annotation file
	#meta-features-  aggregation of a set of features
	#featureCounts program uses the gene_id attribute available in the GTF format
	#General format
		###featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] 	
#Summarize paired-end reads and counting fragments (instead of reads) using a user-provided GTF annotation file:
#Only need to run this step once and then you can analyze the counts matrix multiple ways 
#Example using multiple samples that are paired end
fcounts <-featureCounts(files="*.sorted.bam",isPairedEnd=TRUE,annot.ext="gtffile",isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id")
write.csv(fcounts$counts, file = "/DIRECTORY_PATH/counts.csv")

#Count matrix input= http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input
#reads in the count matrix made in the previous step- whole ocunt matrix
cts <- as.matrix(read.csv("/DIRECTORY_PATH/counts.csv", sep=",",row.names=1,stringsAsFactors=TRUE)

#Alternatively, reads in a subset of the count matrix made in the previous step- only analyzing subset of samples
fcts_path <- file.path(indir, "counts.csv")
fcts <- read.csv(fcts_path)
fcts
library(dplyr)
#List files you will be analyzing- needs to match the nomenclature and order used in the sampletable
fcts_noHS <- fcts %>% select("X", "RAL304_noHS.sorted.bam","RAL335_noHS.sorted.bam","RAL737_noHS.sorted.bam","RAL819_noHS.sorted.bam","RAL195_noHS.sorted.bam","RAL69_noHS.sorted.bam","RAL93_noHS.sorted.bam","RAL359_noHS.sorted.bam","RAL387_noHS.sorted.bam","RAL409_noHS.sorted.bam")
write_csv(fcts_noHS, "/DIRECTORY_PATH/fcts_noHS.csv")
cts <- as.matrix(read.csv("/DIRECTORY_PATH/fcts_noHS.csv", sep=",",row.names=1,stringsAsFactors=TRUE))

#reads in sample information table for DESeq2 analysis
sampleTable_DESeq <- read.csv(sampleTable_DESeq_path, row.names=1)
coldata <- sampleTable_DESeq
	#Example of sampleTable_DESeq.csv
	#,precondGroup
	#RAL304_noHS.sorted.bam,,detrimental
	#RAL335_noHS.sorted.bam,detrimental
	#RAL737_noHS.sorted.bam,detrimental
	#RAL819_noHS.sorted.bam,detrimental
	#RAL195_noHS.sorted.bam,detrimental
	#RAL69_noHS.sorted.bam,beneficial
	#RAL93_noHS.sorted.bam,beneficial
	#RAL359_noHS.sorted.bam,beneficial
	#RAL387_noHS.sorted.bam,beneficial
	#RAL409_noHS.sorted.bam,beneficial

# examine the count matrix and column data to see if they are consistent in terms of sample order
head(cts,2)
coldata
#Absolutely critical that the columns of the count matrix and the rows of the column data are in the same order and use the same naming 
#If they are not in the same order run:
all(rownames(coldata) %in% colnames(cts)) 
all(rownames(coldata) == colnames(cts)) #Needs to say TRUE to continue
ts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts)) #Needs to say TRUE to continue

#Construct a DESeqDataSet
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ precondGroup)
dds
	

##Step 3-Differential expression analysis##

#pre-filtering
#pre-filter low count genes before running the DESeq2 functions- optional but useful
#Removes genes with 0 counts for all samples:
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]

#Factor Levels note
#R will choose a reference level for factors based on alphabetical order
#Note- in this case detrimental is the control 
#Option1: Use Factor 
dds$precondGroup <- factor(dds$precondGroup, levels = c("detrimental","beneficial"))
#Option2-relevel, just specifying the reference level:
dds$precondGroup <- relevel(dds$precondGroup, ref = "detrimental")

#Run differential expression analysis-
dds <- DESeq(dds)
res <- results(dds)
res

#p-values and adjusted p-values
#We can order our results table by the smallest p value
resOrdered <- res[order(res$pvalue),]
summary(res)
#Set FDR cutoff to 10%
res10 <- results(dds, alpha=0.1)
summary(res05)

#Write csv of results
write.csv(as.data.frame(resOrdered), file="condition_treated_results.csv")
#Exporting only the results which pass an adjusted p value threshold
resSig <- subset(resOrdered, padj < 0.1)

##Step4- Data transformations and visualization
#both vst and rlog aim to remove the dependence of the variance on the mean
#Extracting transformed values for plots
#The running times are shorter when using blind=FALSE and if the function DESeq has already been run
	#then it is not necessary to re-estimate the dispersion values.
vsd <- vst(dds, blind=FALSE)

rld <- rlog(dds, blind = FALSE) <- used this data transformation for the final results in the manuscript 
head(assay(rld), 3)
plotPCA(rld, intgroup = c("precondGroup"))
pcaData <- plotPCA(rld, intgroup = c( "precondGroup"), returnData = TRUE)
pcaData


##Plot PCA 

#Perform rlog transformation
rld <- rlog(dds, blind = FALSE)

#Compute PCA using prcomp
pca <- prcomp(assay(rld))

# Plot PCA with group labels:
plot <- plotPCA(rld, intgroup = "precondGroup")

# Customize the plotPCA function to include desired axis limits
library(ggplot2)
pcaData <- plotPCA(rlog, intgroup = precondGroup, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(xlim=c(-10,15), ylim(-10,10))

#Identifying outliers from PCA
#1-Identify outliers visually by observing points that deviate significantly from the main cluster.
    #Outliers are usually those samples that are distant from the majority of the other samples in the plot.
	#You can look for points that are located far away from the main cluster or in a different direction.
#2- To programmatically identify outliers, you can calculate the distances of each sample from the centroid of the PCA plot 
	#and define a threshold to determine outliers. You can use the dist() function to calculate the Euclidean distances between 
	#each sample and the centroid:

pcaData <- plotPCA(rld, intgroup = c( "precondGroup"), returnData = TRUE)

#Matrix of pairwise Euclidean distances between the data points in the PCA space
distances <- dist(pcaData[, c("PC1", "PC2")])
print(distances)

##To calculate the distance of each data point from the center (origin) in the two-dimensional PCA space, you can use the Pythagorean theorem (Euclidean distance formula)
#Calculate the distance of each data point from the center (origin)
distances_from_center <- sqrt(pcaData$PC1^2 + pcaData$PC2^2)
# View the distances
print(distances_from_center)

#Calculate threshold, in this case 2 SD away from the mean 
# Calculate the threshold
threshold <- mean(distances_from_center) + 2 * sd(distances_from_center)
threshold


~~~~
devtools::session_info()
#should be included at end of script to encourage reproducible work. just makes a list of the version of each program you used. Nice to have for your final results. 

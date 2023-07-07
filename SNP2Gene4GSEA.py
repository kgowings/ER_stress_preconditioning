...
#OBJECTIVES
#GWAS_file example (includes header): 
  #chr	rs	ps	n_miss	allele1	allele0	af	beta	se	l_remle	l_mle	p_wald	p_lrt	p_score
  #1	2L_5317_SNP	5317	15	A	G	0.064	-1.147246e-01	4.671726e+00	1.000000e-05	1.000000e+05	9.804403e-01	4.788057e-01	4.802449e-01

##step 1 
#take the log2 of the last column and make it a new column 

#SNP2gene file examples(SNP2geneOUTPUT.txt, no  header)
  #2L_10004508_SNP	SiteClass[FBgn0266167|CR44873|EXON|0;FBgn0043456|CG4747|INTRON|0]
  #2L_10005703_SNP	SiteClass[FBgn0043456|CG4747|SYNONYMOUS_CODING|0;FBgn0266167|CR44873|UPSTREAM|660]
  #2L_10007544_SNP	SiteClass[FBgn0043456|CG4747|INTRON|0]
  #2L_10012262_SNP	SiteClass[FBgn0011584|Trp1|DOWNSTREAM|144;FBgn0032175|CG13131|DOWNSTREAM|965]
  #2L_10024064_SNP	SiteClass[|||]

##Step 2 (Snp2GeneWithPValue_verbose.txt)
  #Make a complied file containing location (second column in GWAS_file, matching first column in SNP2geneOUTPUT file), p_score (from GWAS file), p_score_log2 (take the log2 of the p_score from first file), and all the info between [] in SNP2gene file (SNP2geneOUTPUT file) 
  #Example of output file: (tab delimited)
    #rs p_score p_score_log2 gene
    #2L_5317_SNP 4.802449e-01 4.802449e-01 FBgn0011584|Trp1|DOWNSTREAM|144;FBgn0011584|Trp1|DOWNSTREAM|144;FBgn0032175|CG13131|DOWNSTREAM|965

##Step 3  (Snp2GeneWithPValue_terse.txt)
#If SiteClass[|||], remove that line - no gene associated with SNP
#Extract FBgn:
#If two options extract FBgn using heirarchy  EXON > NON_SYNONYMOUS_CODING > STOP_GAINED > STOP_LOST > CODON_CHANGE_PLUS_CODON_DELETION > SPLICE_SITE_REGION > 
#SYNONYMOUS_CODING > UTR_5_PRIME > UTR_3_PRIME > INTRON > (DOWNSTREAM or UPSTREAM if several than the one with the shortest number following it )

#3: Final output file:
#rs p_score_log2 gene
#2L_10024064_SNP 4.802449e-01 FBgn0011584
'''

import math, decimal

GWAS_file = open('GWAS_file.txt')
GWAS_file.readline() ## trims header
SNP2gene_file = open('SNP2geneOUTPUT.txt')
verbose_output = open('Snp2GeneWithPValue_verbose.txt', 'w')
terse_output = open('Snp2GeneWithPValue_terse.txt', 'w')
terse_output.write('entry_id' + '\t' + 'experiment1' + '\n')
verbose_output.write('location' + '\t' + 'chr' + '\t' + 'rs' + '\t' + 'ps' + '\t' + 'n_miss'+ '\t' + 'allele1' + '\t' + 'allele0' + '\t' + 'af' + '\t' + 'beta' + '\t' + 'se' + '\t' + 'l_remle' + '\t' + 'l_mle' + '\t' + 'p_wald' + '\t' + 'p_lrt' + '\t' + 'p_score' + '\t' + 'p_score_log2' + '\t' + 'geneInfo' + '\n')

print 'Read GWAS and SNP2gene files, starting processing'

GWASByGene = {}
for line in GWAS_file:
    key = line.split('\t')[1]
    if key in GWASByGene:
        print 'Error with this line in GWAS file: ' + line #avoids duplicate entries for the same key 
    else:
        GWASByGene[key] = line # assigns the entire GWAS_file line as the value for that key in the dictionary if it doesn't alread
#Reads GWAS_file and populates a dictionary named 'GWASByGene' using the second tab-separated field as the key and the entire line as the value

print 'Indexed GWAS file. Iterating through gene file'

infoPrioritization = [
    'EXON',
    'NON_SYNONYMOUS_CODING',
    'STOP_GAINED',
    'FRAME_SHIFT',
    'STOP_LOST',
    'START_LOST',
    'CODON_CHANGE_PLUS_CODON_INSERTION',
    'CODON_CHANGE_PLUS_CODON_DELETION',
    'CODON_INSERTION',
    'CODON_DELETION',
    'SPLICE_SITE_REGION',
    'SYNONYMOUS_CODING',
    'SYNONYMOUS_STOP',
    'START_GAINED',
    'NON_SYNONYMOUS_START',
    'UTR_5_PRIME',
    'UTR_3_PRIME',
    'INTRON',
    # DOWNSTREAM and UPSTREAM are always last
]

for line in SNP2gene_file:
    chunks = line.split('\t') #creates a string

    location = chunks[0] #ex 2L_10004508_SNP
    geneInfo = chunks[1].split('SiteClass[')[1].split(']\n')[0] 
      #chunks[1] retrieves the element at index 1 from the chunks list or string, ex SiteClass[FBgn0266167|CR44873|EXON|0;FBgn0043456|CG4747|INTRON|0]
      #split('SiteClass[') -splits the element wherever the 'SiteClass[' occurs and returns a list of substrings
      #[1] to retrieve the second element from the list, ex FBgn0266167|CR44873|EXON|0;FBgn0043456|CG4747|INTRON|0]...
      #split(']\n') - splits new element wherever the substring ']\n' occurs and returns a list of substrings
      #[0] to retrieve the first element from the new list, ex FBgn0266167|CR44873|EXON|0;FBgn0043456|CG4747|INTRON|0
      #extracted substring is assigned to the variable geneInfo
    GWASLine = ''.join(GWASByGene[location].split('\n'))
      #Takes the value associated with the location key in the GWASByGene dictionary 
      #Splits it into a list of substrings using newline characters as delimiters
      #Then joins these substrings together using newline characters as separators to form a single string, ex 1	2L_5317_SNP	5317	15	A	G	0.064	-1.147246e-01	4.671726e+00	1.000000e-05	1.000000e+05	9.804403e-01	4.788057e-01	4.802449e-01
    GWASChunks = GWASLine.split('\t')
    p_score = ''.join(GWASChunks[13].split('\n')) # remove trailing newline
    p_score_log2 = str(-math.log(float(p_score), 2))

    # writing to verbose file
    verbose_output.write(location + '\t' + GWASLine + '\t' + p_score_log2 + '\t' + geneInfo + '\n')

    if ''.join(geneInfo.split('|')) == '': # it's a [|||] line
        continue

    infoChunks = geneInfo.split(';')
    if len(infoChunks) == 0: # no chunks there; skipping
        continue
    terseInfo = None

    for priority in infoPrioritization:
        if terseInfo: # priority already matched
            break
        for chunk in infoChunks:
            try:
                if chunk.split('|')[2] == priority:
                    terseInfo = chunk.split('|')[0]
                    break
            except:
                print 'Failed with infoChunk: ' + chunk
                print '(on line: ' + line + ' )'

    if terseInfo == None: # DOWNSTREAM or UPSTREAM logic
        terseInfo = infoChunks[0].split('|')[0]
        lowestStreamValue = int(infoChunks[0].split('|')[3])
        for chunk in infoChunks:
            if 'DOWNSTREAM' not in chunk and 'UPSTREAM' not in chunk:
                print 'Could not understand this line: ' + line
            val = int(chunk.split('|')[3])
            if val < lowestStreamValue:
                terseInfo = chunk.split('|')[0]

    terse_output.write(terseInfo + '\t' + p_score_log2 + '\n')

print 'Done!'

#Adding header to terse file 
#echo -e "snp_location\tentry_id\texperiment1" | cat - snp2gene_pvalue_terse.txt > snp2gene_pvalue_terse_headers.txt

#grabbing second two columns from terse file 
#  awk 'BEGIN { FS = "\t" } ; { print $2, $3 }' snp2gene_pvalue_terse_headers.txt > snp2gene_pvalue_terse_headers_GSEAready.txt

#convert comma to tab 
#awk '{gsub(/\,/,"\t");print;}' ddv502supp_tables4.txt > GWAS_tab.txt


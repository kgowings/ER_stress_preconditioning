# example line from GWAfile (OUTPUT_FILE_NAME.assoc.txt file):
# 1	2L_19227086_SNP	chr2L:19227086			2	T	G	0.026	4.14E+01	6.96E+00	4.56E+00 ...

# query based on example line: 2L_19227086_SNP

# example line from SNPfile (dgrp.fb557.txt):
# 2L_10000016_SNP	C	SiteClass[FBgn0051755|SoYb|NON_SYNONYMOUS_CODING|0;FBgn0051875|CG31875|INTRON|0],TranscriptAnnot[NON_SYNONYMOUS_CO ...

# section to output from example line:
# SiteClass[FBgn0051755|SoYb|NON_SYNONYMOUS_CODING|0;FBgn0051875|CG31875|INTRON|0]

import sys

GWAfile = open('OUTPUT_FILE_NAME.assoc.txt file')
SNPfile = open('dgrp.fb557.txt')
SNP2gene_OUTPUT = open('SNP2gene_OUTPUT.txt', 'w')

####
#
# create list of queries from GWAfile
#
####

sys.stderr.write('creating list of queries...\n')

# discard header line
GWAfile.readline()

# queries = list()
queries = set()

for line in GWAfile:
    columns = line.split('\t')
    query = columns[1]
    queries.add(query)


####
#
# iterate through SNPfile. if "key" is in the list of queries from GWAfile, output "value"
#
####

sys.stderr.write('searching for queries in SNPfile...\n')

for line in SNPfile:
    try:
        key = line.split('\t')[0]
        value = line.split('\t')[2].split(',')[0]
        if key in queries:
			SNP2gene_OUTPUT.write(key + '\t' + value + '\n')
    except:
        print "problem line:", line
SNP2gene_OUTPUT.close()

sys.stderr.write('Done!\n')


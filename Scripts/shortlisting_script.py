#!/usr/bin/env python

import argparse
import pandas
import statsmodels
from statsmodels.stats.multitest import fdrcorrection
pandas.options.mode.chained_assignment = None

     
# Create the argument parser that takes in a csv file of LD grouping, and the CoDeS3D output files 'eqtls.txt' , 'snps.txt', and 'genes.txt'

parser = argparse.ArgumentParser(description='')
parser.add_argument("-l", "--ld_data", type= argparse.FileType('r'), required=True, metavar="ld_data.txt", help = "A comma-separated file containing the LD groupings of the input SNPs, the first column should contain the risk loci, and the second column the input SNPs", dest = "ld_data")
parser.add_argument("-e","--eqtls", type=argparse.FileType('r'), required=True, metavar = "eqtls.txt", help = "A tab-separated file containing eQTL-target gene pairs and its corresponding P-values, this is the CoDeS3D 'eqtls.txt' output file", dest = "eqtls")
parser.add_argument("-g", "--genes", type=argparse.FileType('r'), required=True, metavar="genes.txt", help = "A tab-separated file showing all the genes that are physically interacting with the input SNPs, this is the CoDeS3D 'genes.txt' output file", dest = "genes")
parser.add_argument("-s", "--snps", type=argparse.FileType('r'), required=True, metavar="snps.txt", help = "A tab-separated file showing all the SNPs are in the eQTL databases, this is the CoDeS3D 'snps.txt' output file", dest = "snps")
args = parser.parse_args()


# Read the input CSV/TSV files into a pandas DataFrame
ld_data = pandas.read_csv(args.ld_data)
eqtls = pandas.read_csv(args.eqtls, sep = '\t')
genes = pandas.read_csv(args.genes, sep = '\t')
snps = pandas.read_csv(args.snps, sep = '\t')

#Create a dictionary that contains the risk loci and its LD partners (risk locus as key, and LD partners as values)
risk_loci_dict = {}

for i in range(len(ld_data)):
    if ld_data.iloc[i,0] not in risk_loci_dict.keys():
        risk_loci_dict[ld_data.iloc[i,0]] = [ld_data.iloc[i,1]]
    else:
        risk_loci_dict[ld_data.iloc[i,0]].append(ld_data.iloc[i,1])



# The function below is to categorize each SNP by the risk locus it belongs to
# For example, "categorize_by_risk_locus('rs2549782')" will output "rs4869313/rs4869314" as the risk locus
def categorize_by_risk_locus(SNP):
    for risk_locus,LD_partners in risk_loci_dict.items():
        if SNP in LD_partners:
            return risk_locus


#At the current state, the "eqtls.txt" file does not contain rsID information for each SNPs, it only contain the variant ID. 
#In order to get the rsID information, we merge "eqtls.txt" with the ['snp','variant_id'] columns of "snps.txt". the ['snp'] column contains the rsID information. 
eqtls =  eqtls.merge(snps[['snp','variant_id']].drop_duplicates(),how='left',left_on='sid',right_on='variant_id')

#At the current state, the "eqtls.txt" file does not contain gene name information for each target gene, it only contain the the gene GENCODE ID. 
#In order to get the gene name information, we merge "eqtls.txt" with the ['gene','gencode_id'] columns of "genes.txt". the ['gene'] column contains the gene name information. 
eqtls = eqtls.merge(genes[['gene','gencode_id']].drop_duplicates(),how='left', left_on='pid', right_on='gencode_id')

#Remove redundant columns from 'eqtls' dataframe
eqtls = eqtls.drop(columns=['sid','pid'])


#Currently, the "eqtls" dataframe already has the SNP rsIDs and target gene name information
#However, the SNPs have not yet been grouped according to the risk locus it belongs to

     
#Now we create a new column called "risk_locus" which will contain the risk locus information of each SNP in the dataframe !
eqtls['risk_locus'] = ''
eqtls['risk_locus'] = eqtls['snp'].map(categorize_by_risk_locus)



#Now, we need to select the spatial eQTL with the lowest p-value for each risk locus-gene-tissue combinations as a representative of the risk locus' regulatory effect on a particular gene in a specific tissue/immune cell type prior to FDR correction.
#The rationale behind this approach is explained in Supplementary Figure 3 in the manuscript
eqtls_minimumP = eqtls.loc[eqtls.groupby(['tissue','risk_locus','gene'])['pval'].idxmin()]

#Save this "eqtls_minimumP" as a csv file called "eqtls_shortlisted.txt". But we want to follow the format of the CoDeS3D output file "eqtls.txt"
eqtls_shortlisted = eqtls_minimumP[['variant_id','gencode_id','sid_chr','sid_pos','adj_pval','pval','b','b_se','maf','tissue']]
eqtls_shortlisted.rename(columns = {'variant_id':'sid','gencode_id':'pid'}, inplace = True)
eqtls_shortlisted.to_csv('eqtls_shortlisted.txt',sep='\t', index=False)


#Now that we have the representative p-value of each hypothesis clusters (i.e., risk locus-gene-tissue combinations)
#We want to do Benjamini-Hochberg FDR correction on the set of representative hypotheses individually in each tissues
#First, we make an empty column called "adj_pval_shortlisted"
eqtls_minimumP['adj_pval_shortlisted'] = ''

#Lets now fill the ['adj_pval_shortlisted'] column with the appropriate value
for x in list(eqtls_minimumP['tissue'].unique()):
    subset = eqtls_minimumP[eqtls_minimumP['tissue'] == x]
    fdr = fdrcorrection(subset['pval'])
    subset['adj_pval_shortlisted']=fdr[1]
    eqtls_minimumP.update(subset)

#Finally.. extract the shortlisted significant eqtls (i.e., rows that have 'adj_pval_shortlisted' less than or equal to 0.05)
significant_eqtls_shortlisted = eqtls_minimumP[eqtls_minimumP['adj_pval_shortlisted'] <= 0.05]


#Re-order the columns for easier reading
significant_eqtls_shortlisted = significant_eqtls_shortlisted.reindex(columns=['risk_locus','variant_id','snp','sid_chr','sid_pos','gencode_id','gene','adj_pval_shortlisted','adj_pval','pval','b','b_se','maf','tissue'])

#Save as a csv file
significant_eqtls_shortlisted.to_csv('significant_eqtls_shortlisted.txt', sep='\t', index=False)






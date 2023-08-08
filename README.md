This repository contains the python scripts used for generating the shortlisted CoDeS3D results for the manuscript "Juvenile idiopathic arthritis-associated genetic loci exhibit spatially constrained gene regulatory effects across multiple tissues and immune cell types". This study combined multiple levels of biological information (3D genome organization [Hi-C] and tissue/cell type-specific eQTL data) to decipher the biological mechanism linking genetic risk to the development of JIA.

The statistical rationale behind this "shortlisting" step is explained below (also in Supplementary Figure 3 in the manuscript)

![Alt text](./shortlisting_illustration.png)

Figure A. Schematic of the method for controlling false discoveries ("shortlisting"). We selected the spatial eQTL with the lowest p-value (indicated with *) in each risk locus-gene- tissue/cell type combination to represent the risk locus' regulatory effect on a particular gene in a specific tissue/immune cell type prior to FDR correction.

In the example outlined in Figure A, left, the lead SNP rs1 (red circle) has six LD SNPs (black circle) making a total of seven SNPs in this rs1-tagged locus. From CoDeS3D, two of these SNPs have potential spatial eQTL interactions to gene A and three SNPs to gene B. This means that we are undertaking five hypothesis tests with five different P-values. However, due to LD structure, the P-values of the neighboring hypotheses are highly correlated, and these hypotheses are likely to be rejected together. At the same time, a rejection of one null hypothesis should not be interpreted as evidence of a causal role for the rejected SNP, rather, we equate discovery with the genomic locus it belongs to. Therefore, we have a discrepancy between the number of rejected hypotheses and the number of discoveries. In this hypothetical scenario, we potentially reject five hypotheses, while making only two biological discoveries (rs1-gene A and rs1-gene B). Due to this mismatch, strategies that control the FDR defined in terms of individual hypotheses (i.e., SNPs) might not lead to satisfactory inference 1. A similar problem has been raised in the literature previously 2–4. Specifically, it has been suggested that in such a situation, a cluster of neighboring correlated rejections should be combined and counted as one rejection. The global false discovery rate should be defined as the "proportion of clusters that are falsely rejected among all rejected clusters" 2. Therefore, here (Supplementary Figure 3, right), we define a cluster of hypotheses in each tissue as all the spatial eQTLs targeting a gene in a particular risk locus. Then, based on the observed data, we selected a representative hypothesis (P-value) of each cluster by selecting the spatial eQTL within the risk locus with the lowest P-value to the indicated target gene, and applying a FDR-controlling procedures to the set of representative hypotheses.

==========================================================================================================================================================================================================================================================================================

In this repository, the "Scripts" folder contains two Python scripts.

1. The first Python script "shortlisting_script.py" takes in three CoDeS3D output files (eqtls.txt, genes.txt, snps.txt), and one csv file of risk locus/LD grouping as the input (see "risk_locus_grouping.csv" file or supplementary table 3 for examples), this script outputs two files :

A). The first output file "significant_eqtls_shortlisted.txt" is a tab-separated file containing spatial eQTL-target gene interactions after shortlisting that passed the FDR correction

B). The second output file "eqtls_shortlisted.txt" is a tab-separated file containing spatial eQTL-target gene interactions after shortlisting that have not been filtered for FDR less than 0.05 

example usage command : "python shortlisting_script.py -e eqtls.txt -g genes.txt -s snps.txt -l risk_locus_grouping.csv"


note :  use "python shortlisting_script.py -h" to access help page








2. However, in the output file "significant_eqtls_shortlisted.txt" , the effect size of each spatial eQTL-target gene interactions is reported with beta value, not allelic fold change (log2[aFc]). In order to get the effect size in term of allelic fold change, we need to use the second python script "summary.py".

the script "summary.py" takes the second output file of the script "shortlisting_script.py" called "eqtls_shortlisted.txt" as the input. Note that in order for "summary.py" to work properly, it must be in the same folder with the CoDeS3D output files "genes.txt" and "snps.txt". The script "summary.py" will output "significant_eqtls.txt" which is a tab-separated file containing spatial eQTL-target gene interactions after shortlisting that passed the FDR correction (effect size is reported using (log2 [aFc]))

example usage command : "python summary.py -e eqtls_shortlisted.txt --eqtl-project GTEx --multi-test tissue -o summary_py_results"


note : the "--multi-test tissue" argument is to make sure that the FDR correction is done separately in each tissues. use "python summary.py -h" to access help page.

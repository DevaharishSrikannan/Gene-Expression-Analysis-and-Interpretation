# Gene-Expression-Analysis-and-Interpretation
Assignment 2 ANAT40040-Bio Principles &amp; Cellular Org
In this study, RNA-seq data is analysed for differential gene expression, survival, and pathway enrichment. Several R programs are used in this investigation for preprocessing, statistical analysis, and visualisation. All of this code's components and their purposes are explained below.

# Step 1 (Download and Extraction)
Initially,Downloaded the TCGA dataset, which includes CNA, clinical, and RNA-seq data, from cBioPortal. To access the appropriate information files, extracted the contents of the tar.gz file into a folder after downloading it. Verify that the extracted files are in the appropriate places for examination.

# Step 2 (Reading the Files and Matching)
To ensure that patient IDs are consistent across the RNA-seq, clinical, and CNA datasets, preprocess the data. By combining gene expression data with clinical and CNA data, this is one of the most significant procedures. To ensure that each and every data item relates to the correct person, combined these datasets according to patient IDs.

# Step 3 (Build the Meta Data)
Build the ERBB2 Amplification Metadata: To Determine whetherthe ERBB2 gene is Amplified or Not Amplified based on the CNA data. It is deemed ERBB2+ if the CNA value for the ERBB2 gene is higher than 0. In order to do additional analysis on the samples, this metadata is necessary for grouping them into ERBB2+ or ERBB2-representative categories.

# Step 4 (Normalization)
DESeq2 corrects the RNA-seq data for biases in library size and sequencing depth. As a result, the results will be similar between samples. A matrix of normalised gene expression data that can be utilised for additional downstream analysis will be the result of this stage.

# Step 5 (Differential Expression Analysis)
 To Identify the overexpression and downregulation status of the genes in the ERBB2+ and ERBB2-groups to conduct differential expression analysis. From the experimental data, log2 fold changes and multiple testing adjustment for p-values was computed. A list of genes with differential expression is the result.

# Step 6 (Pathway Enrichment Analysis)
Pathway Enrichment analysis is the Analysis to test the enriched biological processes or pathways associated with the differentially expressed genes: tools such as Gene Ontology (GO), KEGG and Reactome has been used for the identification of pathways associated with the differentially expressed genes. The Pathway Enrichment analysis aims at revealing some insights into the biological meaning of the differentially expressed genes.

# Step 7 (Variance Stabilized Transformation (VST))
A Variance Stabilizing Transformation is applied to the RNA-seq data to bring down heteroscedasticity and better prepare data for downstream analysis. It stabilizes the variance across genes, making the data more comparable. The transformed data is subject to both PCA and heatmap generation.

# Step 8 (Principal Component Analysis and Heatmap)
PCA and Heatmap Generation: In this method, involves the generation of a PCA plot, which visualizes the global clustering of samples based on their gene expression. Then it is convienent to check if there is any difference between the ERBB2+ and ERBB2- sample groups. Within the generated heatmap, one is able to visualize differential gene expression patterns between the different samples.

# Step 9 (Survival Analysis)
Survival Analysis Using glmnet: The variance-stabilized expression values of differentially expressed genes put into the Cox proportional hazards model. The model showed the effect of gene expression on patient survival because the glmnet package applies regularized regression to this end and, thus, evaluate gene expression effects on survival outcome.

# supervised_vs_unsupervised_comparison

This repository contains the files necessary for R code which will compare the ability of random forest and K-means clustering to distinguish mega populations based on SNP. The mega-populations used are east asian and european. The genes from which the SNP were obtained are the beta subunit of hemoglobin (betaglobin) as well as the pdha1 gene. All VCF and tsv files were obtained from the human genome data slicer (http://grch37.ensembl.org/Homo_sapiens/Tools/DataSlicer?db=core) and the 1000 genomes project (https://www.internationalgenome.org/)

Files contained include:

sup_unsup_compare: The R code that will obtain all SNP data from both mega populations in the VCF files, then classify each individual as belonging to one of the two mega populations based on their SNP patterns using both random forest as well as k-means classification. This code will output a confusion matrix for the random forest classifier to determine its efficacy as well as a silhouette index for the k-means classifier.

betaglob_eas.vcf: A vcf file containing the SNP information for betaglobin from the east asian mega-population
betaglob_eu.vcf: A vcf file containing the SNP information for betaglobin from the european mega-population
pdha1_eas.vcf: A vcf file containing the SNP information for pdha1 from the east asian mega-population
pdha1_eu.vcf: A vcf file containing the SNP information for pdha1 from the european mega-population

igsr_populations.tsv: A tab separated document containing information regarding each mega-population as well as their respective sub-populations

sup_unsup_compare_rmarkdown : A pdf file generated using rmarkdown showing the results that would be expected if the sup_unsup_compare R code was executed.

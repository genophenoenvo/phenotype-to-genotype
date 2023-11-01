<h1>Post-GWAS Prioritization of  Genome-Phenome Associations in Sorghum
</h1>

## Abstract
This study presents an innovative approach for understanding the genetic underpinnings of  two key phenotypes in Sorghum bicolor: maximum canopy height and maximum growth rate. Genome-Wide Association Studies (GWAS) are widely used to decipher the genetic basis of traits in organisms, but the challenge lies in selecting an appropriate statistically significant threshold for analysis. Our goal was to employ GWAS to pinpoint the genetic markers associated with the phenotypes of interest using specific permissive-filtered threshold values that allows the inclusion of broader collections of explanatory candidate genes. Then, we utilized a pattern recognition technique to prioritize a set of informative genes, which hold potential for further investigation and could find applications in Artificial Intelligence systems. Utilizing a subset of the Sorghum Bioenergy Association Panel cultivated at the Maricopa Agricultural Center in Arizona, we sought to unveil patterns between phenotypic similarity and genetic proximity among accessions in order to organize Single Nucleotide Polymorphisms (SNPs) which are likely to be associated with the phenotypic trait. Additionally, we explored the impact of this method by considering all SNPs versus focusing on SNPs classified through the GWAS pre-filter. Experimental results indicated that our approach effectively prioritizes SNPs and genes influencing the phenotype of interest. Moreover, this methodology holds promise for feature selection from genomic data for predicting complex phenotypic traits influenced by numerous genes and environmental conditions, and could pave the way for further research in this field.

## Pipeline of Our Method
![Image not available.](figures/Figure1.jpg)

## Datasets Used:
* _Normalized_ phenotypic trait data for MAC Season 6: https://github.com/genophenoenvo/JAGS-logistic-growth

* Phenotypic trait data (end-of-season height) for Clemson for generalizability test: https://github.com/genophenoenvo/terraref-datasets 

* GWAS-filtered SNPs of _Sorghum Bicolor_ based on phenotypes maximum canopy height and maximum growth rate at various p-values: https://github.com/genophenoenvo/sorghum_data/releases/tag/v0.0.4

* All variants (SNPs) of _Sorghum Bicolor_: https://storage.googleapis.com/gpe-sorghum/whole-vcf-snp-arrays

## Files and Folders:
* **/target_correlation_matrix:** Contains target correlation matrices created based on phenotypic measurements.
  
* **/SNP_outputs_v0.0.4:** Contains the list of candidate and significant SNP IDs extracted based on GWAS-filtered SNPs after execution our algorithm.

* **/SNP_outputs_v0.0.5:** Contains the list of candidate and significant SNP IDs extracted based on all SNPs after execution our algorithm.

* **/significant_SNPs:** TSV expansion of significant SNPs from the folders SNP_outputs_v0.0.4 and SNP_outputs_v0.0.5 showing snpEff annotation information.

* **/candidate_SNPs:** TSV expansion of candidate SNPs from the folders SNP_outputs_v0.0.4 and SNP_outputs_v0.0.5 showing snpEff annotation information.

* **/correlation_trend_v0.0.4:** Containining the plots depicting the correlation trend while adding SNPs to the "candidate" SNP list for the experiments based on GWAS-filtered SNPs.

* **/correlation_trend_v0.0.5:** Containining the plots depicting the correlation trend while adding SNPs to the "candidate" SNP list for the experiments based on all SNPs.

* **/codes:** Contains the codes used in this work.

* **/sorghum_panthar_outputs:**  

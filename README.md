<h1>Phenotype to Genotype: Toward Determining Genetic Basis for Canopy Height Differences in Sorghum </h1>

## Summary
This work describes a novel approach to deduce the genetic basis of two phenotypes of Sorghum bicolor, viz., maximum canopy height and maximum growth rate. A subset of the Sorghum Bioenergy Association Panel was grown at the Maricopa Agricultural Center in Arizona as part of the TERRA-REF project, where various phenotype measurements were taken regularly throughout the growing season by an automated field scanner system and field scientists. In this work, we utilized a logistic growth curve model to normalize the phenotype measurements in order to account for the environmental variations. The differences in normalized phenotypic observations among accessions were used to build a “target" similarity matrix that describes the pairwise proximity of accessions. At the same time, multiple similarity matrices denoting pairwise genetic relatedness of 274 accessions were constructed based on the Single Nucleotide Polymorphisms (SNPs) identified using Genome-wide Association Studies; each of the SNP corresponds to a similarity matrix. Thereafter, we adopted a pattern recognition technique to identify which of these SNP similarity matrices are analogous to the “target” similarity matrix. This approach helped us identify the SNPs that are most likely to impact the phenotype of interest.

## Pipeline of Our Method
![Image not available.](figures/Figure1.jpg)

## Datasets Used:
* _Normalized_ phenotypic trait data for MAC Season 6: https://github.com/genophenoenvo/JAGS-logistic-growth
* Phenotypic trait data (end-of-season height) for Clemson for generalizability test: https://github.com/genophenoenvo/terraref-datasets 
* GWAS-filtered SNPs of _Sorghum Bicolor_ based on phenotypes maximum canopy height and maximum growth rate at varipus p-values: https://github.com/genophenoenvo/sorghum_data/releases/tag/v0.0.4
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

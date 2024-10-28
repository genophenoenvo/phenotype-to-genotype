<h1>Post-GWAS Prioritization of  Genome-Phenome Associations in Sorghum
</h1>

## Abstract
Genome-Wide Association Studies (GWAS) are widely used to infer the genetic basis of traits in organisms, yet selecting appropriate thresholds for analysis remains a significant challenge. In this study, we developed the Sequential SNP Prioritization Algorithm (SSPA) to elucidate the genetic underpinnings of two key phenotypes in Sorghum bicolor: maximum canopy height and maximum growth rate. Utilizing a subset of the Sorghum Bioenergy Association Panel cultivated at the Mar-icopa Agricultural Center in Arizona, our objective was to employ GWAS with specific permis-sive-filtered thresholds to identify the genetic markers associated with these traits, allowing for a broader collection of explanatory candidate genes. Following this, our proposed method incorpo-rates a feature engineering approach based on statistical correlation coefficient to reveal patterns between phenotypic similarity and genetic proximity across 274 accessions. This approach helps prioritize Single Nucleotide Polymorphisms (SNPs) likely to be associated with the studied phe-notype. Additionally, we evaluated the impact of SSPA by considering all variants (SNPs) as inputs, without any GWAS filtering, as a complementary analysis. Empirical evidence including ontolo-gy-based gene function, spatial and temporal expression, and similarity to known homologs, demonstrated that SSPA effectively prioritizes SNPs and genes influencing the phenotype of in-terest, providing valuable insights for functional genetics research.

## Pipeline of Our Method
![Image not available.](figures/Outline.pdf)

## Datasets Used:
* _Normalized_ phenotypic trait data for MAC Season 6: https://github.com/genophenoenvo/JAGS-logistic-growth

* Phenotypic trait data (end-of-season height) for Clemson for generalizability test: https://github.com/genophenoenvo/terraref-datasets 

* GWAS-filtered SNPs of _Sorghum Bicolor_ based on phenotypes maximum canopy height and maximum growth rate at various p-values: https://github.com/genophenoenvo/sorghum_data/releases/tag/v0.0.4

* All variants (SNPs) of _Sorghum Bicolor_: https://storage.googleapis.com/gpe-sorghum/whole-vcf-snp-arrays

## Files and Folders:
* **/target_correlation_matrix:** Contains phenotype similarity matrices created based on phenotypic measurements.
  
* **/SNP_outputs_v0.0.4:** Contains the list of candidate and prioritized SNP IDs extracted based on GWAS-filtered SNPs after execution our algorithm.

* **/SNP_outputs_v0.0.5:** Contains the list of candidate and prioritized SNP IDs extracted based on all SNPs after execution our algorithm.

* **/significant_SNPs:** TSV expansion of prioritized SNPs from the folders SNP_outputs_v0.0.4 and SNP_outputs_v0.0.5 showing snpEff annotation information.

* **/candidate_SNPs:** TSV expansion of candidate SNPs from the folders SNP_outputs_v0.0.4 and SNP_outputs_v0.0.5 showing snpEff annotation information.

* **/correlation_trend_v0.0.4:** Containining the plots depicting the correlation trend while adding SNPs to the "candidate" SNP list for the experiments based on GWAS-filtered SNPs.

* **/correlation_trend_v0.0.5:** Containining the plots depicting the correlation trend while adding SNPs to the "candidate" SNP list for the experiments based on all SNPs.

* **/codes:** Contains the codes used in this work.

* **/sorghum_panthar_outputs:**  

# PUBH-8878---Final-Project
Global Population Structure and Allele Frequency Variation at Breast Cancer Risk Loci


### Environment
- **R version:** 4.5.1  
- **bcftools:** v1.22  
- **bgzip/tabix:** htslib 1.22  

### R Packages Used
| Purpose | Packages |
|----------|-----------|
| Data import and wrangling | `vcfR`, `dplyr`, `tidyr`, `purrr` |
| Visualization | `ggplot2`, `ggrepel` |

You can filter Variant Call Files in your environment using:

bcftools view -r 13:32000000-33000000 \  -Oz -o brca2_region.vcf.gz \ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
bcftools view -r 22:28000000-30000000 \-Oz -o chek2_region.vcf.gz \ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
bcftools view -r 17:7400000-7800000 \-Oz -o tp53_region.vcf.gz \ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

bcftools index brca2_region.vcf.gz
bcftools index chek2_region.vcf.gz
bcftools index tp53_region.vcf.gz

zgrep -v "^#" brca2_region.vcf.gz | head
zgrep -v "^#" chek2_region.vcf.gz | head
zgrep -v "^#" tp53_region.vcf.gz | head

#### Take produced files and run them in R using the PUBH8878.R 

Berdecio, I. (2025). Global Population Structure and Allele Frequency Variation at Breast Cancer Risk Loci. PUBH 8878, George Washington University.

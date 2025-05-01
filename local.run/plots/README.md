## Inputs
ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
## Outputs
This fold store some outputs from our NF pipeline with different parameter setting:
- maf_0.5-rand_0.02.pca_plot.png
  - params.maf_interval = "(0.05, 0.5]" // select VCF records where MAF >0.05
  - params.rand_fraction = 0.02          // select 2% VCF records randomly on top of MAF selection
- rand_0.02.pca_plot.png
  - params.rand_fraction = 0.02          // select 2% VCF records randomly from original vcf (chr22)

- maf_0.3-0.31.pca_plot.png
  - params.maf_interval = "(0.3, 0.31]" // select VCF records where MAF >0.3 and MAF <= 0.31


## Outputs
This fold store some outputs from our NF pipeline. 
- maf_0.5-rand_0.02.pca_plot.png with parameters
  - params.maf_interval = "(0.03, 0.31)" // select VCF records where MAF >0.5
  - params.rand_fraction = 0.02          // select 2% VCF records randomly on top of MAF selection
- rand_0.02.pca_plot.png
  - params.rand_fraction = 0.02          // select 2% VCF records randomly from original vcf (chr22)


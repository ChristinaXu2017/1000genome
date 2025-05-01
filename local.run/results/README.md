
## Experiment with single input
### Inputs
- ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz <br>(Downloaded from s3://1000genomes/release/20130502/)
### Outputs
This fold store some outputs from our NF pipeline with different parameter setting:
- maf_0.5-rand_0.02.pca_plot.png
  - params.maf_interval = "(0.05, 0.5]" <br> ( select VCF records where MAF >0.05 )
  - params.rand_fraction = 0.02   <br> ( select 2% VCF records randomly on top of MAF selection )
- maf_0.3-0.31.pca_plot.png
  - params.maf_interval = "(0.3, 0.31]" <br> ( elect VCF records where MAF >0.3 and MAF <= 0.31 )
- rand_0.02.pca_plot.png
  - params.rand_fraction = 0.02          <br> ( select 2% VCF records randomly from original chr22 vcf )
  - rand_0.02.timeline.html & rand_0.02.trace.txt <br> ( JOb wall time and RAM usage)


## Experiment with multi inputs
### Inputs
- ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz  
- ALL.2000.vcf.gz  ( a testing data with 2000 VCF record only )

### Outputs
- 2inputs.rand.pca_plot.png
  - params.rand_fraction = 0.02          
  - 2inputs.rand.pca_plot.png.timeline.html & 2inputs.rand.pca_plot.png.trace.txt  

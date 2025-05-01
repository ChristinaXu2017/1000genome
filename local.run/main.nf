#!/usr/bin/env nextflow

params.chunk_size = 500
params.reads = "/Users/xu102/Documents/gitHub/xu.repo.nextflow/1000genome/data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
params.panel = "/Users/xu102/Documents/gitHub/xu.repo.nextflow/1000genome/data/integrated_call_samples_v3.20130502.ALL.panel"
params.maf_interval = "(0.03, 0.31)" // select all records, hence skip maf selection
// params.rand_fraction = 0.02


include { SPLIT_SAMPLE } from './modules/scripts'
include { GET_COUNTS } from './modules/scripts'
include { MERGE_COUNTS } from './modules/scripts'
include { CREATE_PCA } from './modules/pca'
include { MAF_SELECT } from './modules/scripts/MAF.nf'
include { RAND_SELECT } from './modules/scripts/MAF.nf'

workflow {


    // Create initial VCF channel 
    in_vcf = params.maf_interval ? 
	MAF_SELECT(Channel.fromPath(params.reads, checkIfExists: true), params.maf_interval).vcf :
 	Channel.fromPath(params.reads, checkIfExists: true)

    if (params.rand_fraction) {
	in_vcf = RAND_SELECT(in_vcf, params.rand_fraction).vcf 
    }

    samples = SPLIT_SAMPLE(in_vcf, params.chunk_size).chunks
    chunk_counts = GET_COUNTS(in_vcf.combine(samples.flatten()))
    all_counts = MERGE_COUNTS(chunk_counts.collect())
    CREATE_PCA(all_counts.counts, channel.fromPath(params.panel))

}



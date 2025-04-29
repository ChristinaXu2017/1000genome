#!/usr/bin/env nextflow

params.chunk_size = 25
params.reads = "/Users/xu102/Documents/gitHub/xu.repo.nextflow/1000genome/data/100sample.2000.vcf"
params.panel = "/Users/xu102/Documents/gitHub/xu.repo.nextflow/1000genome/data/integrated_call_samples_v3.20130502.ALL.panel"
params.maf_interval = "(0, 0.03]" // select all records, hence skip maf selection
params.rand_fraction = 0.1


include { SPLIT_SAMPLE } from './modules/scripts'
include { GET_COUNTS } from './modules/scripts'
include { MERGE_COUNTS } from './modules/scripts'
include { CREATE_PCA } from './modules/pca'
include { MAF_SELECT } from './modules/scripts/MAF.nf'

workflow {


    // Create initial VCF channel 
    in_vcf = params.maf_interval ? 
	MAF_SELECT(Channel.fromPath(params.reads, checkIfExists: true), params.maf_interval).vcf :
 	Channel.fromPath(params.reads, checkIfExists: true)

    samples = SPLIT_SAMPLE(in_vcf, params.chunk_size).chunks
    chunk_counts = GET_COUNTS(in_vcf.combine(samples.flatten()))
    all_counts = MERGE_COUNTS(chunk_counts.collect())
    CREATE_PCA(all_counts.counts, channel.fromPath(params.panel))

}



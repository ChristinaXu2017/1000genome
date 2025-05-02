#!/usr/bin/env nextflow

include { MAF_SELECT } from './modules/scripts/MAF.nf'
include { RAND_SELECT } from './modules/scripts/MAF.nf'
include { SPLIT_SAMPLE } from './modules/scripts'
include { GET_COUNTS } from './modules/scripts'
include { MERGE_COUNTS } from './modules/scripts'
include { CREATE_PCA } from './modules/pca'

workflow {


    // Create initial VCF channel 
    in_vcfs = params.maf_interval ? 
	MAF_SELECT(Channel.fromPath(params.reads, checkIfExists: true), params.maf_interval).vcf :
 	Channel.fromPath(params.reads, checkIfExists: true)

    if (params.rand_fraction) {
	in_vcfs = RAND_SELECT(in_vcfs, params.rand_fraction).vcf 
    }

    // get subsample from one of the input vcf
    samples = SPLIT_SAMPLE(in_vcfs.first(), params.chunk_size).chunks

    chunk_counts = GET_COUNTS(in_vcfs.combine(samples.flatten()))
    all_counts = MERGE_COUNTS(chunk_counts.collect())
    CREATE_PCA(all_counts.counts, channel.fromPath(params.panel))

}



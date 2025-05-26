#!/usr/bin/env nextflow

include { MAF_SELECT } from './modules/bcftool/select.nf'
include { RAND_SELECT } from './modules/bcftool/select.nf'
include { SPLIT_SAMPLE } from './modules/bcftool/count.nf'
include { GET_COUNTS } from './modules/bcftool/count.nf'
include { MERGE_COUNTS } from './modules/bcftool/count.nf'
include { CREATE_PCA } from './modules/pca'

workflow {

    // Create input channel
    in_vcfs = Channel.fromPath(params.reads, checkIfExists: true)

    // Step1: ingesting
    if ( params.rand_fraction ) {
	in_vcfs = RAND_SELECT(in_vcfs, params.rand_fraction).vcf 
    }
    if (params.maf_interval) {
	in_vcfs = MAF_SELECT(in_vcfs, params.maf_interval).vcf 
    }

    // Step2: counting variants parallel. here subsample from one of the input vcf
    samples = SPLIT_SAMPLE(in_vcfs.first(), params.chunk_size).chunks
    chunk_counts = GET_COUNTS(in_vcfs.combine(samples.flatten()))
    all_counts = MERGE_COUNTS(chunk_counts.collect())

    // Step3: plotting
    CREATE_PCA(all_counts.counts, channel.fromPath(params.panel))
}

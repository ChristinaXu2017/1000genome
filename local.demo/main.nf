#!/usr/bin/env nextflow

include { MAF_SELECT } from './modules/bcftool/select.nf'
include { RAND_SELECT } from './modules/bcftool/select.nf'
include { SPLIT_SAMPLE } from './modules/bcftool/count.nf'
include { GET_COUNTS } from './modules/bcftool/count.nf'
include { MERGE_COUNTS } from './modules/bcftool/count.nf'
include { CREATE_PCA } from './modules/pca'

params.chunk_size = 500
params.rand_fraction = 0.02
params.maf_interval = "(0.001, 0.5]"
params.reads = "/Users/xu102/Documents/gitHub/xu.repo.nextflow/1000genome/data/ALL.*.vcf.gz"
params.panel = "/Users/xu102/Documents/gitHub/xu.repo.nextflow/1000genome/data/integrated_call_samples_v3.20130502.ALL.panel"
params.results  = "/Users/xu102/Documents/gitHub/1000genome/local.run/results"


workflow {

    // Create input channel
    in_vcfs = Channel.fromPath(params.reads, checkIfExists: true)

    // Call processes
    if ( params.rand_fraction ) {
	in_vcfs = RAND_SELECT(in_vcfs, params.rand_fraction).vcf 
    }

    if (params.maf_interval) {
	in_vcfs = MAF_SELECT(in_vcfs, params.maf_interval).vcf 
    }

    samples = SPLIT_SAMPLE(in_vcfs.first(), params.chunk_size).chunks

    chunk_counts = GET_COUNTS(in_vcfs.combine(samples.flatten()))
    all_counts = MERGE_COUNTS(chunk_counts.collect())
    CREATE_PCA(all_counts.counts, channel.fromPath(params.panel))
}

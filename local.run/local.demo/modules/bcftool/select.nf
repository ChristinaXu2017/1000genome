#!/usr/bin/env nextflow

process RAND_SELECT {

    container "staphb/bcftools:1.21"
    publishDir "$params.results", mode: 'symlink'

    input:
    	path vcf
    	val select_fraction

    output:
    	path "rand_${vcf}", emit: vcf

    script:
    """
        bcftools view ${vcf} | \
        awk -v frac="${select_fraction}" 'BEGIN {srand()} /^#/ {print} !/^#/ {if (rand() < frac) print}' | \
        bcftools view -o rand_${vcf} -Oz
    """
}

process MAF_SELECT {
    container "staphb/bcftools:1.21"

    input:
        path vcf
        val maf_interval   

    output:
        path "filtered_${vcf}", optional: true, emit: vcf

    script:
    """
        #!/bin/bash

        # Parse lower bound and operator
        if echo "${maf_interval}" | grep -q '\\['; then
            lower_op=">="
            lower_bound=\$(echo "${maf_interval}" | sed 's/.*\\[\\(.*\\),.*/\\1/');
        else
            lower_op=">"
            lower_bound=\$(echo "${maf_interval}" | sed 's/.*(\\(.*\\),.*/\\1/');
        fi

        # Parse upper bound and operator
        if echo "${maf_interval}" | grep -q '\\]'; then
            upper_op="<="
            upper_bound=\$(echo "${maf_interval}" | sed 's/.*,\\s*\\(.*\\)\\].*/\\1/')
        else
            upper_op="<"
            upper_bound=\$(echo "${maf_interval}" | sed 's/.*,\\s*\\(.*\\)).*/\\1/')
        fi

        # Trim leading and trailing spaces
        upper_bound=\$(echo "\$upper_bound" | tr -d '[:space:]')
        lower_bound=\$(echo "\$lower_bound" | tr -d '[:space:]')


        # Validate bounds are numeric
        if ! [[ "\$lower_bound" =~ ^[0-9]*\\.?[0-9]*\$ ]] || ! [[ "\$upper_bound" =~ ^[0-9]*\\.?[0-9]*\$ ]]; then
            echo "Error: Invalid bounds in maf_interval: ${maf_interval}" >&2
            exit 1
        fi

        # Construct bcftools filter expression
        filter_expr="INFO/MAF \${lower_op} \${lower_bound} && INFO/MAF \${upper_op} \${upper_bound}"
        bcftools +fill-tags ${vcf} -Oz -- -t MAF | bcftools view -i "\${filter_expr}" -o filtered_${vcf} -Oz
    """
}



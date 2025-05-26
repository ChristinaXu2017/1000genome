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
	#!/bin/bash

	num_region=100  # bcftool query 100 region 
	contig_id=\$( echo ${vcf} | sed -n 's/.*chr\\([0-9]\\{1,\\}\\).*/\\1/p')
	chr_length=\$(bcftools view -h "${vcf}" | grep "^##contig" | grep "ID=\$contig_id," | awk -F'[= >]' '{print \$5}')

# per region is chr_length/100; window is 2x5% of region
# total 10% of whole chr region will be selected
window=\$((chr_length / num_region * 5 / 100))

# Generate 100 random region using awk
awk -v seed=\$RANDOM -v n=\$num_region -v len=\$chr_length '
BEGIN {
    srand(seed);  # Seed random number generator
    for (i = 1; i <= n; i++) {
        pos = int(rand() * len) + 1;  # Random position from 1 to chrom_length
        print pos;
    }
}' | sort -n | uniq > position.txt

# Select fraction of variants per window
# we have to increase 10 times of fraction, due to window is 10% of entire chr
bcftools index ${vcf}
while read -r start; do
    win_start=\$((start - window))
    win_end=\$((start + window))
    [ \$win_start -lt 0 ] && win_start=0

    bcftools view ${vcf} \${contig_id}:\${win_start}-\${win_end} | \
    awk -v frac="${select_fraction}" 'BEGIN {srand()} /^#/ {print} !/^#/ {if (rand() / 10 < frac) print}' | \
    bcftools view -o tmp.\${start}.${vcf} -Oz
   
    # create index for merge    
    bcftools index tmp.\${start}.${vcf}
done < position.txt

input_files=(tmp.*.${vcf})
if [ \${#input_files[@]} -gt 1 ]; then
    # merge multi tmp vcfs
    bcftools merge --force-samples -o rand_${vcf} tmp.*.${vcf}
else
    # rename the single vcf to final
    mv \${input_files} rand_${vcf}
fi

    """


}


// awk rand() is slow for all variants
process RAND_SELECT1 {

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



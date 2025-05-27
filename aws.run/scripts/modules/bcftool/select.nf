#!/usr/bin/env nextflow

process MAF_SELECT {
    label 'bcftool'

    input:
        path vcf
        // e.g., [0, 0.5) for MAF < 0.5
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


process RAND_SELECT {
    label 'bcftool'

    input:
        path vcf
        val select_fraction

    // emit multi vcf outputs based on region number
    output:
        path "rand.*.${vcf}", emit: vcfs

    script:
    """
        #!/bin/bash

        contig_id=\$( echo ${vcf} | sed -n 's/.*chr\\([0-9]\\{1,\\}\\).*/\\1/p')
        chr_length=\$(bcftools view -h "${vcf}" | grep "^##contig" | grep "ID=\$contig_id," | awk -F'[= >]' '{print \$5}')
        
        # divid chr to regions with size 9m, here chr22 length is 51m
        size_region=9000000  # 9m position is 9/52 of chr22 
        num_region=\$((chr_length / size_region + 1))

        # total 50% of region will be selected
        region_percent=0.5
        window=\$((size_region / 2))
        echo "chr length is \$chr_length; region no is \$num_region; window size is \$window"

        # Generate 100 random region using awk
        awk -v seed=\$RANDOM -v n=\$num_region -v len=\$chr_length '
        BEGIN {
            srand(seed);  # Seed random number generator
            for (i = 1; i <= n; i++) {
                pos = int(rand() * len) + 1;  # Random position from 1 to chrom_length
                print pos;
            }
        }' | sort -n | uniq > position.txt

        # create multi output vcfs with selected fraction of variants per window 
        bcftools index ${vcf}
        last_end=1
        while read -r start; do
            win_start=\$([ \$start -lt \$last_end ] && echo "\$last_end" || echo \$start )
            win_end=\$(( win_start + window ))
            last_end=\$(( win_end + 1 ))
            echo "query \${contig_id}:\${win_start}-\${win_end} ..."

            bcftools view ${vcf} \${contig_id}:\${win_start}-\${win_end} | \
            awk -v frac="${select_fraction}" -v perc="\${region_percent}"  'BEGIN {srand()} /^#/ {print} !/^#/ {if (rand() * perc  < frac) print}' | \
            bcftools view -o rand.\${start}.${vcf} -Oz
        done < position.txt

    """
}


process RAND_SELECT_bf {
    label 'bcftool'
    //container = 'job-definition://arn:aws:batch:us-east-1:041314368896:job-definition/fargate-job-bcftool:2'

    input:
	    path vcf
	    // e.g., 0.3 for 30% of records to select
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


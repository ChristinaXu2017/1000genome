#!/usr/bin/env nextflow

process MERGE_COUNTS {

    publishDir "$params.results/allel_counts", mode: 'copy', overwrite: true

    input:
      path fcount

    output:
      path "allel_counts.csv",  emit: counts

    script:
    """

      fmerge="allel_counts.csv"
      touch \${fmerge}
	
      # Get the first file
      first_file=\$(ls ${fcount} | head -n 1)
      head -n 1 \$first_file > \${fmerge}

      # Append data (skipping header) from all files
      for file in ${fcount}; do
    	tail -n +2 \$file >> \${fmerge}
      done	

    """
}

process GET_COUNTS {

    input:
      tuple path(vcf), path(chunk)  // Expect a tuple of [vcf, chunk]

    output:
      path "${chunk.name}.${vcf.name}.allel_count.csv", emit: chunk

    script:
    """
      output=${chunk.name}.${vcf.name}.allel_count.csv
      # Write CSV header
      echo "sample,variant,counter" > \$output
	
      samples=(); 
      while IFS= read -r line; do samples+=("\$line"); done < ${chunk}
      for sample in "\${samples[@]}"; do
	        bcftools query -s "\$sample" -f '%ID\\t%FORMAT\\n' ${vcf} | grep -v "0|0" | grep -v "0/0"  | \
	        awk -v sample="\$sample" '{
       		    variant = \$1
       		    # Count the number of non-reference alleles
        	    split(\$NF, genotypes, "[|/]")
        	    counter = 0;
        	    for (i in genotypes) {
            		  if (genotypes[i] != "0" && genotypes[i] != "") counter++
        	    }
          	  print sample "," variant "," counter
  	      }' >> \$output

      done
    """
}

process SPLIT_SAMPLE {
     
    // one channel of a file; another channel with a value
    input: 
      path vcf; 
      val no_sample;

    // output a channel named chunks with multi files
    output:
      path "sample.chunk_*", emit: chunks

    script:
    """
      bcftools query -l ${vcf} | split -l ${no_sample} - sample.chunk_
    """
}

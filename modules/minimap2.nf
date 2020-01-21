/*Comment section: */

process minimap2_fasta {
  label 'minimap2'
  publishDir "${params.output}/${name}/", mode: 'copy', pattern: "*.gz" 

  input: 
    tuple val(name), file(fasta)
    file(db)

  output:
    file("*.gz")

  script:
    """
    minimap2 -ax asm5 -t ${task.cpus} -o ${name}.sam ${db} ${fasta}
    samtools fasta -f 4 -0 ${name}.clean.fasta ${name}.sam
    samtools fasta -F 4 -0 ${name}.contamination.fasta ${name}.sam
    gzip -f ${name}.clean.fasta
    gzip -f ${name}.contamination.fasta
    rm ${name}.sam
    """
}

process minimap2_nano {
  label 'minimap2'
  publishDir "${params.output}/${name}/", mode: 'copy', pattern: "*.gz" 

  input: 
    tuple val(name), file(fastq)
    file(db)

  output:
    file("*.gz")

  script:
    """

    # remove spaces in read IDs to keep them in the later cleaned output
    if [[ ${fastq} =~ \\.gz\$ ]]; then
      zcat ${fastq} | sed 's/ /DECONTAMINATE/g' > ${name}.id.fastq
    else
      sed 's/ /DECONTAMINATE/g' ${fastq} > ${name}.id.fastq
    fi

    minimap2 -ax map-ont -t ${task.cpus} -o ${name}.sam ${db} ${name}.id.fastq
    samtools fasta -f 4 -0 ${name}.clean.id.fastq ${name}.sam
    samtools fasta -F 4 -0 ${name}.contamination.id.fastq ${name}.sam

    sed 's/DECONTAMINATE/ /g' ${name}.clean.id.fastq | gzip > ${name}.clean.fastq.gz
    sed 's/DECONTAMINATE/ /g' ${name}.contamination.id.fastq | gzip > ${name}.contamination.fastq.gz
     
    rm ${name}.sam ${name}.clean.id.fastq ${name}.contamination.id.fastq ${name}.id.fastq
    """
}

process minimap2_illumina {
  label 'minimap2'
  publishDir "${params.output}/${name}/ill_martin_extraction", mode: 'copy', pattern: "*.gz" 

  input: 
    tuple val(name), file(r1), file(r2)
    file(db)

  output:
    file("*.gz")

  script:
    """
    # replace the space in the header to retain the full read IDs after mapping (the mapper would split the ID otherwise after the first space) 
    if [[ ${r1} =~ \\.gz\$ ]]; then
      zcat ${r1} | sed 's/ /DECONTAMINATE/g' > ${name}.R1.id.fastq
    else
      sed 's/ /DECONTAMINATE/g' ${r1} > ${name}.R1.id.fastq
    fi
    if [[ ${r2} =~ \\.gz\$ ]]; then
      zcat ${r2} | sed 's/ /DECONTAMINATE/g' > ${name}.R2.id.fastq
    else
      sed 's/ /DECONTAMINATE/g' ${r2} > ${name}.R2.id.fastq
    fi

    # Use samtools -F 2 to discard only reads mapped in proper pair:
    minimap2 -ax sr -t ${task.cpus} -o ${name}.sam ${db} ${r1}.R1.id.fastq ${r2}.R2.id.fastq
    samtools fastq -F 2 -1 ${name}.clean.R1.id.fastq -2 ${name}.clean.R2.id.fastq ${name}.sam
    samtools fastq -f 2 -1 ${name}.contamination.R1.id.fastq -2 ${name}.contamination.R2.id.fastq ${name}.sam

    # restore the original read IDs
    sed 's/DECONTAMINATE/ /g' ${name}.clean.R1.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/1"}else{print \$0};LINE++;}' | gzip > ${name}.clean.R1.fastq.gz 
    sed 's/DECONTAMINATE/ /g' ${name}.clean.R2.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/2"}else{print \$0};LINE++;}' | gzip > ${name}.clean.R2.fastq.gz
    sed 's/DECONTAMINATE/ /g' ${name}.contamination.R1.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/1"}else{print \$0};LINE++;}' | gzip > ${name}.contamination.R1.fastq.gz 
    sed 's/DECONTAMINATE/ /g' ${name}.contamination.R2.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/2"}else{print \$0};LINE++;}' | gzip > ${name}.contamination.R2.fastq.gz

    # remove intermediate files
    rm ${name}.R1.id.fastq ${name}.R2.id.fastq ${name}.clean.R1.id.fastq ${name}.clean.R2.id.fastq ${name}.contamination.R1.id.fastq ${name}.contamination.R2.id.fastq ${name}.sam

    """
}

process minimap2_illumina_ebi_extraction {
  label 'minimap2'
  publishDir "${params.output}/${name}/ill_ebi_extraction", mode: 'copy', pattern: "*.gz" 

  input: 
    tuple val(name), file(r1), file(r2)
    file(db)

  output:
    file("*.gz")

  script:
    """
    # replace the space in the header to retain the full read IDs after mapping (the mapper would split the ID otherwise after the first space) 
    if [[ ${r1} =~ \\.gz\$ ]]; then
      zcat ${r1} | sed 's/ /DECONTAMINATE/g' > ${name}.R1.id.fastq
    else
      sed 's/ /DECONTAMINATE/g' ${r1} > ${name}.R1.id.fastq
    fi
    if [[ ${r2} =~ \\.gz\$ ]]; then
      zcat ${r2} | sed 's/ /DECONTAMINATE/g' > ${name}.R2.id.fastq
    else
      sed 's/ /DECONTAMINATE/g' ${r2} > ${name}.R2.id.fastq
    fi

    # Use samtools -F 2 to discard only reads mapped in proper pair:
    minimap2 -ax sr -t ${task.cpus} -o ${name}.sam ${db} ${r1}.R1.id.fastq ${r2}.R2.id.fastq
    samtools fastq -f 12 -F 256 -1 ${name}.clean.R1.id.fastq -2 ${name}.clean.R2.id.fastq ${name}.sam
    samtools fastq -f 2 -1 ${name}.contamination.R1.id.fastq -2 ${name}.contamination.R2.id.fastq ${name}.sam

    # restore the original read IDs
    sed 's/DECONTAMINATE/ /g' ${name}.clean.R1.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/1"}else{print \$0};LINE++;}' | gzip > ${name}.clean.R1.fastq.gz 
    sed 's/DECONTAMINATE/ /g' ${name}.clean.R2.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/2"}else{print \$0};LINE++;}' | gzip > ${name}.clean.R2.fastq.gz
    sed 's/DECONTAMINATE/ /g' ${name}.contamination.R1.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/1"}else{print \$0};LINE++;}' | gzip > ${name}.contamination.R1.fastq.gz 
    sed 's/DECONTAMINATE/ /g' ${name}.contamination.R2.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/2"}else{print \$0};LINE++;}' | gzip > ${name}.contamination.R2.fastq.gz

    # remove intermediate files
    rm ${name}.R1.id.fastq ${name}.R2.id.fastq ${name}.clean.R1.id.fastq ${name}.clean.R2.id.fastq ${name}.contamination.R1.id.fastq ${name}.contamination.R2.id.fastq ${name}.sam

    """
}
process filter_un_mapped_alignments {
  label 'minimap2'

  input:
    tuple val(name), path(sam), path(reads)
    val(mode)

  output:
    tuple val(name), val('clean'), path('*clean.fast{q,a}'), emit: cleaned_reads
    tuple val(name), val('contamination'), path('*contamination.fast{q,a}'), emit: contaminated_reads

  script:
  if ( mode == 'paired' ) {
    """
    # Use samtools -F 2 to discard only reads mapped in proper pair:
    samtools fastq -F 2 -1 ${reads[0].baseName}.clean.fastq -2 ${reads[1].baseName}.clean.fastq ${name}.sam
    samtools fastq -f 2 -1 ${reads[0].baseName}.contamination.fastq -2 ${reads[1].baseName}.contamination.fastq ${name}.sam
    """
  } else if ( mode == 'single' || mode == 'fasta' ) {
    dtype = (mode == 'single') ? 'q' : 'a'
    """
    samtools fast${dtype} -f 4 -0 ${reads.baseName}.clean.fast${dtype} ${sam}
    samtools fast${dtype} -F 4 -0 ${reads.baseName}.contamination.fast${dtype} ${sam}
    """
  } else {
    error "Invalid mode: ${mode}"
  }
}

process make_contamination_bam {
  label 'minimap2'

  publishDir "${params.output}/${name}/${tool}", mode: 'copy', pattern: "*.contamination.sorted.bam*"

  input:
    tuple val(name), path(sam), path(reads)
    val(mode)
    val(tool)

  output:
    tuple val(name), path ('*.contamination.sorted.bam'), emit: contamination_bam
    tuple val(name), path ('*.contamination.sorted.bam.bai'), emit: contamination_bai
    tuple val(name), path ('idxstats.tsv'), emit: idxstats

  script:
  if ( mode == 'paired' ) {
    """
    samtools view -b -f 2 -F 2048 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
    samtools index ${name}.contamination.sorted.bam
    samtools idxstats ${name}.contamination.sorted.bam > idxstats.tsv
    """
  } else {
    """
    samtools view -b -F 2052 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
    samtools index ${name}.contamination.sorted.bam
    samtools idxstats  ${name}.contamination.sorted.bam > idxstats.tsv
    """
  }
}

process filter_soft_clipped_alignments {
  label 'samclipy'
  label 'smallTask'

  publishDir "${params.output}/${name}/${tool}", mode: 'copy', pattern: "*.bam*"

  input:
  tuple val(name), path (bam)
  val (minClip)
  val (tool)
  
  output:
  path ('*.ambiguous.bam')
  path ('*.ambiguous.bam.bai')
  
  script:
  """
  git clone https://github.com/MarieLataretu/samclipy.git --branch v0.0.2 || git clone git@github.com:MarieLataretu/samclipy.git --branch v0.0.2 
  samtools view -h ${bam} | python samclipy/samclipy.py --invert --minClip ${minClip} | samtools sort > ${name}.ambiguous.bam
  samtools view -h ${bam} | python samclipy/samclipy.py --minClip ${minClip} | samtools sort > ${name}.unambiguous.bam
  samtools index ${name}.ambiguous.bam
  samtools index ${name}.unambiguous.bam
  """
}

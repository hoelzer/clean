/*Comment section: */

process minimap2_fasta {
  label 'minimap2'

  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.contamination.sorted.bam*"
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*clean.fasta.gz"
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*contamination.fasta.gz"

  input: 
    tuple val(name), path(fasta)
    path db

  output:
    tuple val(name), env(TOTALCONTIGS), emit: num_contigs
    tuple val(name), path('*.sam'), path(fasta), emit: sam // reads just for naming

  script:
  """
  if [[ ${fasta} =~ \\.gz\$ ]]; then
    TOTALCONTIGS=\$(zgrep '^>' ${fasta} | wc -l)
  else
    TOTALCONTIGS=\$(grep '^>' ${fasta} | wc -l)
  fi

  minimap2 -ax asm5 -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${fasta}s
  """
}

process minimap2_nano {
  label 'minimap2'

  input: 
    tuple val(name), path(reads)
    path (db)

  output:
    tuple val(name), path('*.sam'), path(reads), emit: sam // reads just for naming

  script:
  """
  PARAMS="-ax map-ont"
  if [[ ${params.reads_rna} != 'false' ]]; then
    PARAMS="-ax splice -k14"
  fi

  minimap2 \$PARAMS -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${reads}
  """
}

process minimap2_illumina {
  label 'minimap2'

  input: 
    tuple val(name), path(reads)
    path (db)
    val (mode)

  output:
    tuple val(name), path('*.sam'), path(reads), emit: sam // reads just for naming

  script:
  if ( mode == 'paired' ) {
    """
    minimap2 -ax sr -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${reads[0]} ${reads[1]}
    """
  } else {
    """
    minimap2 -ax sr -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${reads}
    """
  }
}

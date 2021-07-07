process compress_reads {
  label 'basics'

  publishDir "${params.output}/${name}/${tool}", mode: 'copy', pattern: "*.gz"

  input:
  tuple val(name), val(type), path(reads)
  val(mode)
  val(tool)

  output:
  tuple val(name), val(type), path("${name}*.${type}.fast{q,a}.gz")

  script:
  if ( mode == 'paired' ) {
    """
    pigz -fc -p ${task.cpus} ${reads[0]} > ${name}_1.${type}.fastq.gz 
    pigz -fc -p ${task.cpus} ${reads[1]} > ${name}_2.${type}.fastq.gz
    """
  } else if ( mode == 'single' || mode == 'fasta' ) {
    dtype = (mode == 'single') ? 'q' : 'a'
    """
    pigz -fc -p ${task.cpus} ${reads} > ${name}.${type}.fast${dtype}.gz
    """
  } else {
    error "Invalid mode: ${mode}"
  }
}

process get_number_of_reads {
  label 'smallTask'

  input:
  tuple val(name), path(reads)
  val(mode)

  output:
  tuple val(name), env(TOTALREADS), emit: totalreads

  script:
  if ( mode == 'paired' ) {
    """
    if [[ ${reads[0]} =~ \\.gz\$ ]]; then
      TOTALREADS_1=\$(zcat ${reads[0]} | echo \$((`wc -l`/4)))
      TOTALREADS_2=\$(zcat ${reads[1]} | echo \$((`wc -l`/4)))
    else
      TOTALREADS_1=\$(cat ${reads[0]} | echo \$((`wc -l`/4)))
      TOTALREADS_2=\$(cat ${reads[1]} | echo \$((`wc -l`/4)))
    fi
    TOTALREADS=\$(( TOTALREADS_1+TOTALREADS_2 ))
    """
  } else {
    """
    if [[ ${reads} =~ \\.gz\$ ]]; then
      TOTALREADS=\$(zcat ${reads} | echo \$((`wc -l`/4)))
    else
      TOTALREADS=\$(cat ${reads} | echo \$((`wc -l`/4)))
    fi
    """
  }
}

process minimap2Stats {
  label 'smallTask'
  
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "stats.txt" 

  input:
  tuple val(name), path(idxstats), val (totalreads)

  output:
  tuple val(name), path ('stats.txt')
  path("${name}_minimap2_stats.tsv"), emit: tsv

  script:
  """
  MAPPEDSUM=\$(awk -F '\\t' '{sum += \$3} END {print sum}' idxstats.tsv)
  UNPROMAPPEDSUM=\$(awk -F '\\t' '/^[^*]/ {sum += \$4} END {print sum}' idxstats.tsv)
  PROPMAP=\$((\$MAPPEDSUM-\$UNPROMAPPEDSUM))

  MAP=\$(awk -v map=\$MAPPEDSUM -v unpromap=\$UNPROMAPPEDSUM -v tot=${totalreads} 'BEGIN {perc=(map-unpromap)/tot*100; print map-unpromap " ("  perc " %) reads were properly mapped; of these:"}')

  FA=\$(awk -v tot=${totalreads} -F '\\t' '/^[^*]/ {propmap=\$3-\$4; if (propmap != 0) print "\\t\\t" propmap " (" propmap/tot*100  "%) reads aligned to " \$1}' idxstats.tsv;)

  touch stats.txt
  cat <<EOF >> stats.txt
  ${totalreads} reads in total; of these:
  \t\$MAP
  \$FA
  EOF

  touch ${name}_minimap2_stats.tsv
  cat <<EOF >> ${name}_minimap2_stats.tsv
  Sample Name\tTotal reads\tMapped reads
  ${name}\t${totalreads}\t\$PROPMAP
  EOF
  """
}

process bbdukStats {
  label 'smallTask'

  publishDir "${params.output}/${name}/bbduk", mode: 'copy', pattern: "stats.txt"

  input:
  tuple val(name), path (bbdukStats)

  output:
  tuple val(name), path ('stats.txt')
  path("${name}_bbduk_stats.tsv"), emit: tsv

  script:
  """
  TOTAL=\$(grep '#Total' ${bbdukStats} | awk -F '\\t' '{print \$2}')
  MNUM=\$(grep '#Matched' ${bbdukStats} | awk -F '\\t' '{print \$2}')
  MPER=\$(grep '#Matched' ${bbdukStats} | awk -F '\\t' '{print \$3}')

  FA=\$(awk -F '\\t' '/^[^#]/ {print "\\t\\t"\$2" ("\$3") aligned to "\$1}' ${bbdukStats})

  touch stats.txt
  cat <<EOF >> stats.txt
  \$TOTAL reads in total; of these:
  \t\$MNUM (\$MPER) reads were properly mapped; of these:
  \$FA
  EOF

  touch ${name}_bbduk_stats.tsv
  cat <<EOF >> ${name}_bbduk_stats.tsv
  Sample Name\tTotal reads\tMapped reads
  ${name}\t\$TOTAL\t\$MNUM
  EOF
  """
}

process writeLog {
  label 'smallTask'

  publishDir "${params.output}/${name}/${tool}", mode: 'copy', pattern: "log.txt"
  
  input:
    val name
    val tool
    val reads
    val db

  output:
    path 'log.txt'
  
  script:
  """
  touch log.txt
  cat <<EOF >> log.txt
  Input reads:\t${reads}
  Contamination:\t${db}
  
  Statistics summary:\t${params.output}/${name}/${tool}/stats.txt
  EOF
  """

}

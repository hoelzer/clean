process download_host {
  label 'minimap2'

  if (params.cloudProcess) {
    publishDir "${params.databases}/hosts", mode: 'copy', pattern: "*.fa.gz" 
  }
  else {
    storeDir "${params.databases}/hosts"
  }

  input:
  val host

  output:
  path "${host}.fa.gz"

  script:
  """
  if [ $host == 'hsa' ]; then
    wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    zcat *.gz | bgzip -@ ${task.cpus} -c > ${host}.fa.gz
  fi
  if [ $host == 'mmu' ]; then
    wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
    zcat *.gz | bgzip -@ ${task.cpus} -c > ${host}.fa.gz
  fi
  if [ $host == 'cli' ]; then
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/337/935/GCF_000337935.1_Cliv_1.0/GCF_000337935.1_Cliv_1.0_genomic.fna.gz
    zcat *.gz | bgzip -@ ${task.cpus} -c > ${host}.fa.gz
  fi
  if [ $host == 'csa' ]; then
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/409/795/GCF_000409795.2_Chlorocebus_sabeus_1.1/GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic.fna.gz
    zcat *.gz | bgzip -@ ${task.cpus} -c > ${host}.fa.gz
  fi
  if [ $host == 'gga' ]; then
    wget ftp://ftp.ensembl.org/pub/release-99/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz
    zcat *.gz | bgzip -@ ${task.cpus} -c > ${host}.fa.gz
  fi
  if [ $host == 'eco' ]; then
    wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//fasta/bacteria_90_collection/escherichia_coli_k_12/dna/Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa.gz
    zcat *.gz | bgzip -@ ${task.cpus} -c > ${host}.fa.gz
  fi
  """
}

process check_own {
  label 'minimap2'

  input:
  path fasta

  output:
  path '*.gz', includeInputs: true

  script:
  """
  # -L for following a symbolic link
  if ! ( file -L $fasta | grep -q 'gzip compressed' ); then
    bgzip -@ ${task.cpus} < ${fasta} > ${fasta}.gz
    # now $fasta'.gz'
  else
    mv ${fasta} ${fasta}.tmp
    zcat ${fasta}.tmp | bgzip -@ ${task.cpus} -c > ${fasta}.gz
  fi
  """
}

process concat_contamination {
  label 'minimap2' // sets conda env (and cpus)
  label 'smallTask' // overrides cpus
  
  publishDir "${params.output}/${name}/${tool}", mode: 'copy', pattern: "db.fa.gz"
  publishDir "${params.output}/${name}/${tool}", mode: 'copy', pattern: "db.fa.fai"

  input:
  val name
  val tool
  path '*'

  output:
  path 'db.fa.gz', emit: fa
  path 'db.fa.fai', emit: fai
  
  script:
  """
  cat * > db.fa.gz
  samtools faidx db.fa.gz
  mv db.fa.gz.fai db.fa.fai
  """
}
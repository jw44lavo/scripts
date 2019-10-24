nextflow.preview.dsl=2

params.mode = "bad" //options: "bad", "verybad" or "nice"
params.path_genomes = "/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/ref"

params.coverage = "30"
params.length = "10000,8500"

params.outdir_badreads = "/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/badreads/${params.mode}/${params.coverage}x_coverage"

/*
* construction of bad reads for nanopore assembly from ref genome with Badread-tool
*
* settings: --reference   ref.fasta
*           --quantitiy   INTEGERx
*           --read length mean,stdev
*           for further settings see (https://github.com/rrwick/Badread#requirements)
*/

/* 
* Fragen:
*         - wie kann ich direkt auf samples zu greifen? mein workflow funktioniert zwar,
            aber ist nicht automatisiert, was ja der sinn dahinter sein sollte...
          - ich nehme aktuell für neue prozesse die files aus meinem wunsch ordner, 
            aber ich könnte sie doch schneller irgendwie über channel aus dem work ordner bekommen oder?
*         - brauche ich überhaupt einen output ordner oder speichere ich das dann doppelt im "work verzeichnis"?
*         - wie kann ich meine OneDrive wirklich aktuell halten? änderungen der datein bermerkt onedrive nicht... 
*           muss ich dafür dann doch einfach git nehmen?
*         - kann ich das work verzeichnis von nextflow verändern?
*/

process create_bad_reads {
  publishDir "${params.outdir_badreads}", mode: "move"
  
  input:
    val path
    each sample

  output:
    file "${sample}_badreads_bad.fastq.gz"

  script:
    if( params.mode == "bad")
      """
      badread simulate --reference ${path}/${sample}.fasta --quantity ${params.coverage}x --length ${params.length} \
      | gzip > ${sample}_badreads_bad.fastq.gz
      """ 


    else if( mode == "verybad")
      """
      badread simulate --reference ${path}/${sample}.fasta --quantity ${params.coverage}x --glitches 1000,100,100 \
      --junk_reads 5 --random_reads 5 --chimeras 10 --identity 75,90,8 --length ${params.length} \
      | gzip > ${sample}_badreads_verybad.fastq.gz
      """

    else if( mode == "nice")
      """
      badread simulate --reference ${path}/${sample}.fasta --quantity ${params.coverage}x --error_model random \
      --qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 --length ${params.length} \
      --chimeras 0 --identity 95,100,4 --start_adapter_seq "" --end_adapter_seq "" \
      | gzip > ${sample}_badreads_nice.fastq.gz
      """   
}


process gunzip {
  publishDir "${params.outdir_badreads}", mode: "copy"
  input:
    val path
    val sample
  
  script:
    """
    gunzip -k ${path}/${sample}_badreads_${params.mode}.fastq.gz
    """
}


process concatenate {
  publishDir "${params.outdir_badreads}", mode: "copy"

  input:
    val path
    val sample
  
  output:
    file "zymo_badreads_${params.mode}_${params.coverage}x.fastq"
  
  script:
  """
  cat ${path}/${sample}_badreads_${params.mode}.fastq >> zymo_badreads_${params.mode}_${params.coverage}x.fastq
  """
}


workflow {
  main:
    samples = Channel.from("c_neoformans") //"l_fermentum, "s_aureus", "l_monocytogenes", "e_faecalis", "b_subtilis", "e_coli", "s_enterica", "p_aeruginosa", "s_cerevisiae", 
    create_bad_reads(params.path_genomes, samples)
    //gunzip(params.outdir_badreads, samples)
    //concatenate(params.outdir_badreads, samples) konkatinieren mit >> überschreibt den vorherigen inhalt

  }
 
/*
  my samples:             badreads 30x:      concatenated 30x:
    l_fermentum               y               n
    s_aureus                  y               n
    l_monocytogenes           y               n
    e_faecalis                y               n
    b_subtilis                y               n
    e_coli                    y               n
    s_enterica                y               n
    p_aeruginosa              y               n
    s_cerevisiae              r               n
    c_neoformans              r               n

    "l_fermentum", "s_aureus", "l_monocytogenes", "e_faecalis", "b_subtilis", "e_coli", "s_enterica", "p_aeruginosa", "s_cerevisiae", "c_neoformans"
*/
nextflow.preview.dsl=2

params.mode = "bad" //options: "bad", "verybad" or "nice"

// params.path_genomes = "/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/ref"

params.coverage = "0.1"
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


/********************************************************************************
* Aufruf: nextflow run create_bad_reads.nf -w /home/johann/Local/work --mode=bad
*********************************************************************************/


/*
basename --> name des files
gunzip(create_..(input_channel))
br.out.collect für alles gleichzeitig
*/
 

process create_bad_reads {

  input:
    file ref_fasta
  
  output:
    file "${ref_fasta.baseName}"
  
  script:
    if (params.mode =="bad")
      """
      badread simulate --reference ${ref_fasta} --quantity ${params.coverage}x --length ${params.length} \
      > ${ref_fasta.baseName}
      """
    
    else if (params.mode == "verybad")
      """
      badread simulate --reference ${ref_fasta} --quantity ${params.coverage}x --glitches 1000,100,100 \
      --junk_reads 5 --random_reads 5 --chimeras 10 --identity 75,90,8 --length ${params.length} \
      > ${ref_fasta.baseName} 
      """

    else if (params.mode == "nice")
      """
      badread simulate --reference ${ref_fasta} --quantity ${params.coverage}x --error_model random \
      --qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 --length ${params.length} \
      --chimeras 0 --identity 95,100,4 --start_adapter_seq "" --end_adapter_seq "" \
      > ${ref_fasta.baseName}
      """   
}


process concatenate {
  publishDir "${params.outdir_badreads}", mode: "copy"

  input:
    file x
  
  output:
    file "zymo_badreads_${params.mode}_${params.coverage}x.fastq"
  
  script:
    """
    cat ${x} >> zymo_badreads_${params.mode}_${params.coverage}x.fastq
    """
}

workflow {
  main:
    params.input = "/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/ref/*.fasta"
    reference_dataset = Channel.fromPath(params.input)
    
    concatenate(create_bad_reads(reference_dataset))

  }
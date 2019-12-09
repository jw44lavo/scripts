nextflow.preview.dsl=2

dir_scripts = "/home/johann/OneDrive/Projects/commonplace/scripts"
dir_env = "/home/johann/OneDrive/Projects/commonplace/environments"



process get_few_reads {
  //publishDir "${dir_scripts}", mode: "copy"

  input:
    file x //all reads
    val u  //number_of_reads
    val v  //current_read_count
  
  output:
    file "${x.baseName}_${u}_${v}.fastq"
  
  script:
    """
    python3 ${dir_scripts}/get_few_reads.py ${x} ${u} ${v}
    """
}

process get_single_reads {
  //publishDir "${dir_scripts}", mode: "copy"

  input:
    file x //input read file

  output:
    file "${x.baseName}_*.fastq"

  script:
    """
    python3 ${dir_scripts}/get_single_reads.py ${x}
    """
}

process get_ids {
  //publishDir "${dir_scripts}", mode: "copy"

  input:
    file x //reference fasta file

  output:
    file "${x.baseName}_ids.txt"

  script:
    """
    grep ">" ${x} | cut -c2- > ${x.baseName}_ids.txt
    """
}

process minimap2_indexing {
  //publishDir "${dir_scripts}", mode: "copy"

  input:
    file x

  output:
    file "${x.baseName}.mmi"

  script:
    """
    minimap2 -d ${x.baseName}.mmi ${x}
    """
}

process minimap2_mapping {
  //publishDir "${dir_scripts}", mode: "copy"

  input:
    file x  //indexed reference genome 
    each y  //read to map

  output:
    file "${y.baseName}_minimap2.sam"

  script:
    """
    minimap2 -ax map-ont ${x} ${y} > ${y.baseName}_minimap2.sam
    """
}

process minimap2_evaluate {
  publishDir "${dir_scripts}", mode: "copy"

  input:
    file x

  output:
    file "${x.baseName}_stats.out"

  script:
    """
    samtools flagstat ${x} >> ${x.baseName}_stats.out
    """
}

process mashmap_mapping {
  //publishDir "${dir_scripts}", mode: "copy"
  
  input:
    file x //reference genome
    each y //read to map

  output:
    file "${y.baseName}_mashmap.out"
  
  script: //Achtung: bei --pi 75 ist der pc eingefroren
    """
    mashmap -r ${x} -q ${y} -s 2500 --pi 85 -o ${y.baseName}_mashmap.out 
    """
}

process mashmap_evaluate {
  publishDir "${dir_scripts}", mode: "copy"

  input:
    file x //mashmap.out file
    file y //reads file --> nur für read count zur berechnung der statistik. !!!noch mal überdenken, da das zählen bei großen files extrem lang dauern könnte!
    file z //id file

  output:
    file "${x.baseName}_stats.out"

  script:
    """
    python3 ${dir_scripts}/evaluate_mapping_mashmap.py ${x} ${y} ${z} >> ${x.baseName}_stats.out
    """
}

process centrifuge_prepare_taxonomy {
  
  conda "${dir_env}/centrifuge2019-11-29.yml"
   
  output:
    file "taxonomy"
  
  script:
    """
    centrifuge-download -o taxonomy taxonomy
    """
}

process centrifuge_prepare_library {
  
  conda "${dir_env}/centrifuge2019-11-29.yml"
  
  output:
    file "seqid2taxid.map"
  
  script:
    """
    centrifuge-download -o library -m -d "bacteria" refseq > seqid2taxid.map
    """
    //-m --> mask low-complexity regions
}

process centrifuge_prepare_reference {
  
  conda "${dir_env}/centrifuge2019-11-29.yml"
  
  input:
    path x

  output:
    file "input-sequences.fna"
  
  script:
    """
    cat ${x}/*/*.fna > input-sequences.fna
    """
}

process centrifuge_indexing {
  
  conda "${dir_env}/centrifuge2019-11-29.yml"
  
  input:
    file seqid2taxid.map 
    file nodes.dmp
    file names.dmp
    file input-sequences.fna

  output:
    file "abv"
  
  script:
    """
    centrifuge-build -p 4 --conversion-table ${seqid2taxid.map} \
                 --taxonomy-tree ${nodes.dmp} --name-table ${names.dmp} \
                 ${input-sequences.fna} abv
    """
    //-p 4 --> 4 Threads
}



process centrifuge_mapping {
  //publishDir "${dir_scripts}", mode: "copy" 
  
  conda "${dir_env}/centrifuge2019-11-29.yml"

  input:
    file x
    each y 

  output:
    stdout "${y.baseName}.out"

  script:
    """
    centrifuge -f -x test ${y}
    """
}



workflow {
    main:

      ref         = Channel.fromPath("/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/genomes/ref/zymo_all/zymo_all.fasta")
      reads       = Channel.fromPath("/home/johann/Local/data/Zymo-GridION-EVEN-BB-SN-PCR-R10HC-flipflop.fq")
      taxonomy_cf = Channel.fromPath("/home/johann/OneDrive/Projects/commonplace/scripts/taxonomy")
      library_cf  = Channel.fromPath("/home/johann/OneDrive/Projects/commonplace/scripts/library")
            
      //r = get_single_reads(few_reads).flatten()

      /********** Data acquisition - Workflow ****************************************************/ 
      few_reads       = get_few_reads(reads, 100, 0)
      ids             = get_ids(ref)


      /********** Minimap2 - Workflow ************************************************************/ 
      
      
      mm2_db          = minimap2_indexing(ref)   //db = Channel.fromPath("${dir_scripts}/zymo_all.mmi")
      mm2_mapping     = minimap2_mapping(mm2_db, few_reads)
      mm2_evaluation  = minimap2_evaluate(mm2_mapping)
      

      /********** MashMap - Workflow *************************************************************/
      
      mmp_mapping     = mashmap_mapping(ref, few_reads)
      mmp_evaluation  = mashmap_evaluate(mmp_mapping, few_reads, ids)


      /********** Centrifuge - Workflow **********************************************************/
      
      //centrifuge_prepare_taxonomy()   funktioniert, dauert kurz (1-3 min)
      //tax_cf = centrifuge_prepare_library()    noch am Testen, dauert sehr lange (~16 h überschlagen)
      //cat_cf = centrifuge_prepare_reference(library_cf)
      //index_cf = centrifuge_indexing(tax_cf, "${taxonomy}/nodes.dmp", "${taxonomy}/names.dmp", cat_cf)
      //centrifuge_mapping(index_cf, r)

}
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
    file "${y.baseName}.sam"

  script:
    """
    minimap2 -ax map-ont ${x} ${y} > ${y.baseName}.sam
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

      ref = Channel.fromPath("/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/genomes/ref/zymo_all/zymo_all.fasta")
      reads = Channel.fromPath("/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/genomes/reads/zymo-GridION_few_reads.fastq")
      taxonomy_cf = Channel.fromPath("/home/johann/OneDrive/Projects/commonplace/scripts/taxonomy")
      library_cf = Channel.fromPath("/home/johann/OneDrive/Projects/commonplace/scripts/library")
            
      r = get_single_reads(reads).flatten()



      // Minimap2 - Workflow
      /*
      db = minimap2_indexing(ref)
      //db = Channel.fromPath("${dir_scripts}/zymo_all.mmi")
      
      minimap2_mapping(db, r)
      */



      // MashMap - Workflow
      /*
      mashmap_mapping(ref, r)
      //mashmap_mapping.out.view()
      */



      // Centrifuge - Workflow
      
      //centrifuge_prepare_taxonomy()   funktioniert, dauert kurz (1-3 min)
      //tax_cf = centrifuge_prepare_library()    noch am Testen, dauert sehr lange (~16 h überschlagen)
      //cat_cf = centrifuge_prepare_reference(library_cf)
      //index_cf = centrifuge_indexing(tax_cf, "${taxonomy}/nodes.dmp", "${taxonomy}/names.dmp", cat_cf)
      //centrifuge_mapping(index_cf, r)

}











/*
process make_diamond_db {

  input:
    file x  //protein database
  
  output:
    file "${x.baseName}.dmnd" //indexed protein database
  
  script:
    """
    diamond makedb --in ${x} --db "${x.baseName}.dmnd"
    """
}

process map_diamond {
  //publishDir "/home/johann/OneDrive/Projects/commonplace/scripts/", mode: "copy"

  input:
   file x //indexed protein database
   file y //reads to map

  output:
    file "${y}_matches.m8" //output in BLAST tabular format
  
  script:
    """
    diamond blastx --db ${x} -q ${y} -o ${y}_matches.m8
    """
}
*/


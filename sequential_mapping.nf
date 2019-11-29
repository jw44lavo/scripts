nextflow.preview.dsl=2

dir_scripts = "/home/johann/OneDrive/Projects/commonplace/scripts"
dir_env = "/home/johann/OneDrive/Projects/commonplace/environments"



process get_few_reads {

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
  publishDir "${dir_scripts}", mode: "copy"
  
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

workflow {
    main:

      ref = Channel.fromPath("/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/genomes/ref/zymo_all/zymo_all.fasta")
      reads = Channel.fromPath("/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/genomes/reads/zymo-GridION_few_reads.fastq")
      
      r = get_single_reads(reads).flatten()
      //r = Channel.fromPath("${dir_scripts}/*.fastq")
     
      /*
      db = minimap2_indexing(ref)
      //db = Channel.fromPath("${dir_scripts}/zymo_all.mmi")
      
      minimap2_mapping(db, r)
      */

      mashmap_mapping(ref, r)
      //mashmap_mapping.out.view()
    
      

}


/*
Idee:
    reads von sequenzierung gruppieren (get_few_reads.py) z.B. 10 (auch mit 1 oder 100 möglich)
    diese n ausgewählten reads mappen
    mapping analysieren und taxonomischen baum herunterwandern
        wiederhole solange mit nächsten gruppen von reads bis bestimmte
        parameter erreicht, sodass aussage über taxon/organismus/stamm möglich
        bzw.
        in "real time" anzeigen, welche entscheidung akutell getroffen

Anmerkungen zur Idee:
    (Pan-)Proteom-Datenbank schon fertig und schon indexed
        (je nach mapper gucken, welche vorverarbeitung nötig)

    externer "read-count" nötig, um nicht immer wieder die gleichen reads zu mappen

*/

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


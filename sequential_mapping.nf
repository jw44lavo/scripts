nextflow.preview.dsl=2

dir_proteom = "/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/proteomes/uniprot/ref/zymo_all"
db_proteom_fasta = "${dir_proteom}/zymo_all_ref_panproteom.fasta"
//db_proteom_indexed = "${dir_proteom}/zymo_all_ref_panproteom.dmnd"

input_reads_loman = "/home/johann/Local/data/Zymo-GridION-EVEN-BB-SN-PCR-R10HC-flipflop.fq"
input_reads_badread = "/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/genomes/badreads/bad/10x_coverage/zymo_badreads_bad_10x.fastq"
params.number_of_reads = 3
current_read_count = 2
reads_format_loman = "fastq" // "fastq" or "fasta"

dir_scripts = "/home/johann/OneDrive/Projects/commonplace/scripts"
dir_env = "/home/johann/OneDrive/Projects/commonplace/environments"


println(" \n \n Welcome to the ultimate sepsis-fast-detection-hero - 'SepSeq' \n \n ")

process get_few_reads {
  //publishDir "${dir_scripts}", mode: "copy"
  input:
    file x //all reads
    val u  //number_of_reads
    val v  //current_read_count
  output:
    file "y.${reads_format_loman}"
  script:
    """
    python3 ${dir_scripts}/get_few_reads.py ${x} y.${reads_format_loman} ${u} ${v}
    """
}

process indexing_diamond {
  publishDir "${dir_scripts}", mode: "copy"
  conda "${dir_env}/mapper2019-11-13.yml"
  input:
    //file x //protein database
  
  output:
    //file "${x.baseName}.dmnd" //indexed protein database 
  
  script:
    """
    diamond makedb --in /home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/proteomes/uniprot/ref/zymo_all/zymo_all_ref_panproteom.fasta -d /home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/proteomes/uniprot/ref/zymo_all/zymo_all_ref_panproteom.dmnd
    """
}

process mapping_diamond {
  input:
   file x //indexed protein database
   file y //reads to map

  output:
    file z //output in BLAST tabular format
  
  script:
  """
  diamond blastx -d nr -q reads.fna -o matches.m8
  """
}

workflow {
    main:
      //get_few_reads(Channel.fromPath(input_reads_loman), params.number_of_reads, current_read_count)
      //indexing_diamond(Channel.fromPath("/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/proteomes/uniprot/ref/zymo_all/zymo_all_ref_panproteom.fasta"))
      //mapping_diamond(indexing_diamond(db_proteom_fasta), get_few_reads(input_reads_loman, params.number_of_reads, current_read_count))

      indexing_diamond() // warum geht das nicht!!? (unten stehender aufruf funktioniert direkt im Terminal)
      /*

      diamond makedb --in /home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/proteomes/uniprot/ref/zymo_all/zymo_all_ref_panproteom.fasta 
      -d /home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/proteomes/uniprot/ref/zymo_all/zymo_all_ref_panproteom.dmnd

      */
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

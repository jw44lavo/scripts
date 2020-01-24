nextflow.preview.dsl=2

//include "./nextflow.config" params(params)

//lib: which holds all the process { ... } specifications

dir_scripts    = "path to scripts"


process get_few_reads {
  publishDir "${dir_scripts}", mode: "copy"

  input:
    file x //all reads
    val u  //wanted number of reads
    val v  //where to start to count
  
  output:
    file "${x.baseName}_${u}_${v}.fastq"
  
  script:
    """
    python3 ${dir_scripts}/get_few_reads.py ${x} ${u} ${v}
    """
}

process get_single_reads {
  publishDir "${dir_scripts}", mode: "copy"

  input:
    file x //input read file

  output:
    file "${x.baseName}_*.fastq"

  script:
    """
    python3 ${dir_scripts}/get_single_reads.py ${x}
    """
}

process count_reads {
  publishDir "${dir_scripts}", mode: "copy"
  
  input:
    file x //reads file

  output:
    file "${x.baseName}_count.out"

  script:
    """
    wc -l ${x} > ${x.baseName}_count.out
    """
}

process get_ids {
  publishDir "${dir_scripts}", mode: "copy"

  input:
    file x //reference fasta file

  output:
    file "${x.baseName}_ids.out"

  script:
    """
    grep "^>" ${x} | cut -c2- > ${x.baseName}_ids.out
    """
}

process minimap2_indexing {
  publishDir "${dir_scripts}", mode: "copy"

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
  publishDir "${dir_scripts}", mode: "copy"

  input:
    file x  //indexed reference genome 
    each y  //reads to map

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
    file x //minimap2.sam file
    file y //reads file --> nur für read count zur berechnung der statistik. !!!noch mal überdenken, da das zählen bei großen files extrem lang dauern könnte!
    file z //id file

  output:
    file "${x.baseName}_stats.out"

  script:
    """
    python3 ${dir_scripts}/evaluate_mapping_minimap2.py ${x} ${y} ${z} > ${x.baseName}_stats.out
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
    python3 ${dir_scripts}/evaluate_mapping_mashmap.py ${x} ${y} ${z} > ${x.baseName}_stats.out
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
    file x

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
    file x 
    file y
    file z
    file i

  output:
    file "abv"
  
  script:
    """
    centrifuge-build -p 4 --conversion-table ${x} \
                 --taxonomy-tree ${y} --name-table ${z} \
                 ${i} abv
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

process mapping_evaluate {
  publishDir "${dir_scripts}", mode: "copy"

  input:
    file x //mapping out file
    file y //reads file --> nur für read count zur berechnung der statistik. !!!noch mal überdenken, da das zählen bei großen files extrem lang dauern könnte!
    file z //id file
    val f  //mapping format
  output:
    file "${x.baseName}_stats.out"

  script:
    """
    python3 ${dir_scripts}/evaluate_mapping.py -i ${x} -c ${y} -n ${z} -o ${x.baseName}_stats.out -f ${f}
    """
}

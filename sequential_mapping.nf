nextflow.preview.dsl=2

dir_scripts = "/home/johann/OneDrive/Projects/commonplace/scripts"
dir_env = "/home/johann/OneDrive/Projects/commonplace/environments"
params.ref = "/home/johann/OneDrive/Projects/commonplace/data/zymo_mock_community/genomes/ref/zymo_all/zymo_all.fasta"
params.reads = "/home/johann/Local/data/Zymo-GridION-EVEN-BB-SN-PCR-R10HC-flipflop.fq"
params.taxonomy_cf = "/home/johann/OneDrive/Projects/commonplace/scripts/taxonomy"
params.library_cf = "/home/johann/OneDrive/Projects/commonplace/scripts/library"
params.mode = "" //minimap2, mashmap, centrifuge, data (for only data acquisition)
params.amount = 100
params.conscious = "no" //yes, no

            

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

process mapping_evaluate {
  publishDir "${dir_scripts}", mode: "copy"

  input:
    file x //mapping out file
    file y //reads file --> nur für read count zur berechnung der statistik. !!!noch mal überdenken, da das zählen bei großen files extrem lang dauern könnte!
    file z //id file

  output:
    file "${x.baseName}_stats.out"

  script:
    """
    python3 ${dir_scripts}/evaluate_mapping.py ${x} ${y} ${z} > ${x.baseName}_stats.out
    """
}

log.info "##############################################################################################"

workflow {
    main:

      ref         = Channel.fromPath(params.ref)
      reads       = Channel.fromPath(params.reads)
      taxonomy_cf = Channel.fromPath(params.taxonomy_cf)
      library_cf  = Channel.fromPath(params.library_cf)
            
      //r = get_single_reads(few_reads).flatten()

      if(params.amount >= 10000 && params.conscious == "no"){
        log.info "uff... thats a big amount of data! are you sure, you want to run the pipeline anyway?"
        log.info "if yes: use '--conscious yes'"
        log.info "if not: see you next time with less data"
      }


      // TRY NESTED SWITCHES!!!

      check_amount = false

      if(params.amount >= 20){
        check_amount = true
      }
      switch(check_amount) {
        case true:
          log.info "uff... thats a big amount of data! are you sure, you want to run the pipeline anyway?"
          log.info "if yes: use '--conscious yes'"
          log.info "if not: see you next time with less data"
          break

        case false:
          log.info "low enough"
          break
        
        defaul:
          break
      }

      switch(params.mode){
        case "":
          log.info "no mapping mode given! please provide '--mode MODE' (minimap2, mashmap, centrifuge,data)"
          break
      
        case "data":
          few_reads       = get_few_reads(reads, params.amount, 0)
          ids             = get_ids(ref)
          read_count      = count_reads(few_reads)
          break

        case "minimap2":
          few_reads       = get_few_reads(reads, params.amount, 0)
          ids             = get_ids(ref)
          read_count      = count_reads(few_reads)
          mm2_db          = minimap2_indexing(ref)
          //mm2_db          = Channel.fromPath("${dir_scripts}/zymo_all.mmi")
          mm2_mapping     = minimap2_mapping(mm2_db, few_reads)
          //mm2_mapping     = Channel.fromPath("${dir_scripts}/Zymo-GridION-EVEN-BB-SN-PCR-R10HC-flipflop_10_0_minimap2.sam")
          //mm2_evaluation  = mapping_evaluate(mm2_mapping, read_count, ids)
          break

        case "mashmap":
          few_reads       = get_few_reads(reads, params.amount, 0)
          ids             = get_ids(ref)
          read_count      = count_reads(few_reads)
          mmp_mapping     = mashmap_mapping(ref, few_reads)
          //mmp_evaluation  = mapping_evaluate(mmp_mapping, read_count, ids)
          break

        case "centrifuge":
          log.info "centrifuge workflow not implemented yet"
          //few_reads       = get_few_reads(reads, params.amount, 0)
          //ids             = get_ids(ref)
          //read_count      = count_reads(few_reads)
          //centrifuge_prepare_taxonomy()   funktioniert, dauert kurz (1-3 min)
          //tax_cf = centrifuge_prepare_library()    noch am Testen, dauert sehr lange (~16 h überschlagen)
          //cat_cf = centrifuge_prepare_reference(library_cf)
          //index_cf = centrifuge_indexing(tax_cf, "${taxonomy}/nodes.dmp", "${taxonomy}/names.dmp", cat_cf)
          //centrifuge_mapping(index_cf, r)
          break

        default:
          log.info "unknown mode! please provide '--mode MODE' (minimap2, mashmap, centrifuge, data)"
      }
}
log.info "##############################################################################################"

nextflow.preview.dsl=2

include "./lib.nf"
include "./nextflow.config" params(params)

//main: exclusively specifies the workflow { ... } and imports params/ processes etc.


log.info "###############################################################################"

workflow {
    main:

      ref         = Channel.fromPath(params.ref)
      reads       = Channel.fromPath(params.reads)
      tax_cf      = Channel.fromPath(params.taxonomy_cf)
      lib_cf      = Channel.fromPath(params.library_cf)
            
      //r = get_single_reads(few_reads).flatten()

      check_amount = false
      if(params.amount <= 10000 || params.amount > 10000 && params.conscious == "yes"){
        check_amount = true
      }

      switch(check_amount){
        case false:
          log.info "uff... thats a big amount of data! are you sure, you want to run the pipeline anyway?"
          log.info "if yes: use '--conscious yes'"
          log.info "if not: see you next time with less data"
          break
      
        case true:
          log.info "let's work!"
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
              mm2_mapping     = minimap2_mapping(mm2_db, few_reads)
              mm2_evaluation  = mapping_evaluate(mm2_mapping, read_count, ids, "mm2")
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
              /*
              few_reads       = get_few_reads(reads, params.amount, 0)
              ids             = get_ids(ref)
              read_count      = count_reads(few_reads)
              //centrifuge_prepare_taxonomy()             //sollte vor der Analyse vorbereitet sein
              //tax_cf = centrifuge_prepare_library()     //sollte vor der Analyse vorbereitet sein
              cat_cf = centrifuge_prepare_reference(library_cf)
              index_cf = centrifuge_indexing("${dir_scripts}/seqid2taxid", "${taxonomy_cf}/nodes.dmp", "${taxonomy_cf}/names.dmp", cat_cf)
              //centrifuge_mapping(index_cf, r)
              */
              break

            default:
              log.info "unknown mode! please provide '--mode MODE' (minimap2, mashmap, centrifuge, data)"
          }
        break
      }

}
log.info "###############################################################################"

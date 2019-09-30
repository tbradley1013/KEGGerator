# function to create keggtap object


keggtap <- function(pathway, match_strict = FALSE, kegg_enzyme = NULL, kegg_orthology = NULL, kegg_pathway = NULL){
  if (is.null(kegg_enzyme)){
    kegg_enzyme <- KEGGerator::kegg_enzymes
  } else if (!is_kegg_tbl(kegg_enzyme, "enzyme")){
    stop("kegg_enzyme must be a kegg_tbl with columns enzyme and enzyme_id", call. = FALSE)
  }

  if (is.null(kegg_orthology)){
    kegg_orthology <- KEGGerator::kegg_orthologies
  } else{
    if (!is_kegg_tbl(kegg_orthology, "orthology")){
      stop("kegg_orthology must be a kegg_tbl with columns orthology and orthology_id", call. = FALSE)
    }
  }

  if (is.null(kegg_pathway)){
    kegg_pathway <- KEGGerator::kegg_pathways
  } else {
    if (!is_kegg_tbl(kegg_pathway, "pathway")){
      stop()
    }
  }

  if (!match_strict){
    pathway_match <- paste(tolower(pathway), collapse = "|")
    pathways <- KEGGerator::kegg_pathways[stringr::str_detect(KEGGerator::kegg_pathways$pathway, pathway_match), ]
  } else {
    pathway_match <- tolower(pathway)
    pathways <- KEGGerator::kegg_pathways[KEGGerator::kegg_pathways]
  }


}

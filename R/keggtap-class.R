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
      stop("kegg_pathway must be a kegg_tbl with columns pathway and pathway_id", call. = FALSE)
    }
  }

  if (!match_strict){
    pathway_match <- paste(tolower(pathway), collapse = "|")
    pathways <- kegg_pathway[stringr::str_detect(tolower(kegg_pathway$pathway), pathway_match), ]
  } else {
    pathway_match <- tolower(pathway)
    pathways <- kegg_pathway[tolower(kegg_pathway$pathway) %in% pathway_match, ]
  }

  if (nrow(pathways) == 0){
    stop("There are no pathways that match your search", call. = FALSE)
  }

  pathway_enzymes <- link_paths(pathways$pathway_id, "enzyme", kegg_enzyme)
  pathway_orthology <- link_paths(pathways$pathway_id, "orthology", kegg_orthology)

  out <- structure(
    list(
      pathway = pathway,
      enzyme = pathway_enzymes,
      orthology = pathway_orthology
    )
  )


}


link_paths <- function(path_id, type, kegg_tbl){

  out <- purrr::map_dfr(path_id, ~{
    path <- stringr::str_replace(.x, "path:", "")

    links <- kegg_link_safe(type, path)

    if (type == "enzyme"){
      out <- tibble::tibble(
        pathway_id = names(links),
        enzyme_id = links
      )
    } else if (type == "orthology"){
      out <- tibble::tibble(
        pathway_id = names(links),
        orthology_id = links
      )
    } else if (type == "module"){
      out <- tibble::tibble(
        pathway_id = names(links),
        module_id = links
      )
    }

    return(out)
  })


  if (type == "enzyme"){
    out <- out %>%
      dplyr::left_join(kegg_tbl, by = "enzyme_id")
  } else if (type == "orthology"){
    out <- out %>%
      dplyr::left_join(kegg_tbl, by = "orthology_id")
  } else if (type == "module"){
    out <- out %>%
      dplyr::left_join(kegg_tbl, by = "module_id")
  }

  return(out)
}


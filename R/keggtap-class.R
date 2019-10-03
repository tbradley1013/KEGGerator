# function to create keggtap object


keggtap <- function(pathway_name, match_strict = FALSE, kegg_enzyme = NULL,
                    kegg_orthology = NULL, kegg_module = NULL,
                    kegg_pathway = NULL, strict = FALSE){

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

  pathway_enzymes <- get_pathway_enzymes(pathway_name, kegg_enzyme = kegg_enzyme, kegg_pathway = kegg_pathway, strict = strict)
  pathway_orthology <- get_pathway_orthologies(pathway_name, kegg_orthology = kegg_orthology, kegg_pathway = kegg_pathway, strict = strict)
  pathway_modules <- get_pathway_modules(pathway_name, kegg_module = kegg_module, kegg_pathway = kegg_pathway, strict = strict)

  out <- structure(
    list(
      pathway = pathway,
      enzyme = pathway_enzymes,
      orthology = pathway_orthology,
      module = pathway_modules
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


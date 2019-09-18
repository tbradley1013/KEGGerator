#' Filter OTUs and kegg organisms
#'
#'
filter_orgs <- function(data, uncertainty, pathway_name, pathways){
  UseMethod("filter_orgs")
}


filter_orgs_internal <- function(org_ids, orgs_tbl, uncert_tbl, uncertainty = 1, pathway_name = NULL, pathways = NULL){

  out <- filter_orgs_uncert(orgs_id = orgs_id, orgs_tbl = orgs_tbl,
                            uncert_tbl = uncert_tbl, uncertainty = uncertainty)

  if (!is.null(pathway_name)){
    out <- filter_orgs_pathway(orgs_id = out, pathway_name = pathway_name, pathways = pathways)
  }

  return(out)

}

filter_orgs_uncert <- function(orgs_id, orgs_tbl, uncert_tbl, uncertainty = 1){

  if (!is_orgs_id(orgs_id)) stop("orgs_id must be of class orgs_id", call. = FALSE)
  if (!is_uncert_tbl(uncert_tbl)) stop("uncert_tbl must be of class uncert_tbl", call. = FALSE)

  otu_keep <- uncert_tbl %>%
    dplyr::filter(total_uncert <= uncertainty) %>%
    dplyr::pull(otu_id)

  otu_keep <- orgs_tbl %>%
    dplyr::filter(otu_id %in% otu_keep)

  org_keep <- otu_keep %>%
    dplyr::pull(genome)

  id_keep <- orgs_id %>%
    dplyr::filter(genome %in% org_keep)

  output <- list(otu_keep = otu_keep, id_keep = id_keep)

  return(output)


}



filter_orgs_pathway <- function(orgs_id, pathway_name, pathways = NULL){
  if (is.null(pathways)){
    pathways <- kegg_pathways
  } else{
    if (!is_kegg_tbl(pathways, "pathway")){
      stop("pathways must be a kegg_tbl with the columns `pathway` and `pathway_id`", call. = FALSE)
    }
  }

  if (!any(stringr::str_detect(pathways$pathway, pathway_name))) stop("pathway_name had no matches in pathways dataset", call. = FALSE)

  pathway_ids <- pathways %>%
    dplyr::filter(stringr::str_detect(pathway, pathway_name)) %>%
    dplyr::mutate(pathway_id = stringr::str_replace(pathway_id, "path:map", "")) %>%
    dplyr::pull(pathway_id)

  output <- orgs_id %>%
    dplyr::mutate(
      genome_id = stringr::str_replace(genome_id, "gn:", ""),
      pathway = purrr::map(genome_id, ~{
        kegg_link_safe("pathway", .x)
      }),
      # pathway = purrr::map(pathway, ~purrr::discard(.x, function(x){!stringr::str_detect(x, pathway_ids)}))
      pathway = purrr::map(pathway, ~{
        prefix <- stringr::str_replace_all(.x, "[0-9]", "") %>%
          unique()

        paths <- stringr::str_replace_all(.x, "[A-z]|\\:", "") %>%
          .[. %in% pathway_ids] %>%
          unique()

        output <- paste0(prefix, paths)
        return(output)
      })
    ) %>%
    dplyr::filter(purrr::map_lgl(pathway, ~{length(.x) > 0}))

  return(output)

}






#' Filter out genomes by pathway
#'
#' @param data a tibble that has genome_ids - likely the output from get_genome_id()
#' @param pathway_name the name of the pathway that you wish to filter
#' @param pathways a kegg_tbl dataset of available kegg pathways. This dataset
#' is produced from the get_kegg_pathway function. If it is NULL (default) the
#' KEGGerator::kegg_pathways dataset will be used
#'
#' @export
filter_genome_by_pathway <- function(data, pathway_name, pathways = NULL){
  if (is.null(pathways)){
    pathways <- kegg_pathways
  } else{
    if (!is_kegg_tbl(pathways, "pathway")){
      stop("pathways must be a kegg_tbl with the columns `pathway` and `pathway_id`", call. = FALSE)
    }
  }

  pathway_ids <- pathways %>%
    dplyr::filter(stringr::str_detect(pathway, pathway_name)) %>%
    dplyr::mutate(pathway_id = stringr::str_replace(pathway_id, "path:map", "")) %>%
    dplyr::pull(pathway_id)

  output <- data %>%
    dplyr::mutate(
      genome_id = stringr::str_replace(genome_id, "gn:", ""),
      pathway = purrr::map(genome_id, ~{
        kegg_link_safe("pathway", .x)
      }),
      # pathway = purrr::map(pathway, ~purrr::discard(.x, function(x){!stringr::str_detect(x, pathway_ids)}))
      pathway = purrr::map(pathway, ~{
        prefix <- stringr::str_replace_all(.x, "[0-9]", "") %>%
          unique()

        paths <- stringr::str_replace_all(.x, "[A-z]|\\:", "") %>%
          .[. %in% pathway_ids] %>%
          unique()

        output <- paste0(prefix, paths)
        return(output)
      })
    ) %>%
    dplyr::filter(purrr::map_lgl(pathway, ~{length(.x) > 0}))

  return(output)
}

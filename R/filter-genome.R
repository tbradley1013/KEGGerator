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
      pathway = map(pathway, ~{
        prefix <- str_replace_all(.x, "[0-9]", "") %>%
          unique()

        paths <- str_replace_all(.x, "[A-z]|\\:", "") %>%
          .[. %in% pathway_ids] %>%
          unique()

        output <- paste0(prefix, paths)
        return(output)
      })
    ) %>%
    dplyr::filter(purrr::map_lgl(pathway, ~{length(.x) > 0}))

  return(output)
}

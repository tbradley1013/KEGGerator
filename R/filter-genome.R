#' Filter out genomes by pathway
#'
#' @param data a tibble that has genome_ids - likely the output from get_genome_id()
#'
#' @export
filter_genome_by_pathway <- function(data, pathway_name){
  pathway_ids <- kegg_pathways %>%
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

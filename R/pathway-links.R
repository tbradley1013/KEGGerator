#' Query all modules associated with a given pathway
#'
#' @param pathway_name a character string specifying a certain pathway. Regular
#' expressions are used to match the pathway so full pathway names do not have
#' to be given
#'
#' @details the value passed to pathway_name is searched against the pathway
#' column in the kegg_pathways dataset
#'
#' @export
get_pathway_modules <- function(pathway_name) {
  pathway <- kegg_pathways %>%
    dplyr::filter(stringr::str_detect(pathway, pathway_name)) %>%
    dplyr::mutate(
      module = purrr::map(pathway_id, ~{
        mods <- KEGGREST::keggLink("module", .x)

        output <- tibble::tibble(
          module_id = mods
        ) %>%
          dplyr::left_join(kegg_modules, by = "module_id")

        return(output)
      })
    ) %>%
    tidyr::unnest(module)

  return(pathway)
}

## THIS IS USED IN THE get_genome_enzymes AND get_genome_orthologies FUNCTIONS
## AND WILL NEED TO BE MODIFIED ONCE CLASSES ARE ADDED TO THESE OUTPUTS IF THAT
## ROUTE IS PURSUED


#' Query all enzymes associated with a given pathway
#'
#' @param pathway_name a character string specifying a certain pathway. Regular
#' expressions are used to match the pathway so full pathway names do not have
#' to be given.
#'
#' @details the value passed to pathway_name is searched against the pathway
#' column in the kegg_pathways dataset
#'
#' @export
get_pathway_enzymes <- function(pathway_name) {
  pathway <- kegg_pathways %>%
    dplyr::filter(stringr::str_detect(pathway, pathway_name)) %>%
    dplyr::mutate(
      enzyme = purrr::map(pathway_id, ~{
        enzymes <- KEGGREST::keggLink("enzyme", .x)

        output <- tibble::tibble(
          enzyme_id = enzymes
        ) %>%
          dplyr::left_join(kegg_enzymes, by = "enzyme_id")

        return(output)
      })
    ) %>%
    tidyr::unnest(enzyme)

  return(pathway)
}

#' Query all orthologies associated with a given pathway
#'
#' @param pathway_name a character string specifying a certain pathway. Regular
#' expressions are used to match the pathway so full pathway names do not have
#' to be given.
#'
#' @details the value passed to pathway_name is searched against the pathway
#' column in the kegg_pathways dataset
#'
#' @export
get_pathway_orthologies <- function(pathway_name){
  pathway <- kegg_pathways %>%
    dplyr::filter(stringr::str_detect(pathway, pathway_name)) %>%
    dplyr::mutate(
      orthology = purrr::map(pathway_id, ~{
        orth <- KEGGREST::keggLink("orthology", .x)

        output <- tibble::tibble(
          orthology_id = orth
        ) %>%
          dplyr::left_join(kegg_orthologies, by = "orthology_id")

        return(output)
      })
    ) %>%
    tidyr::unnest(orthology)

  return(pathway)
}

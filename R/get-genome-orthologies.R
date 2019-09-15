#' Query all orthologies related to each genome in the dataset
#'
#' @param data a tibble with genomes and genome_id or the output from `get_genome_id()`
#' @param pathway_orthologies a vector of pathway orthologies
#' @param kegg_orthology a kegg_tbl with the columns orthology and orthology_id
#' which can be generated with the get_kegg_orthology function.
#'
#' @details if `pathway_orthologies` is NULL (default) than all orthologies for each genome
#' will be returned. If a vector of orthology_ids is provided than only the orthologies related
#' to each genome that are in this vector will be returned. This is useful if used in
#' coordination with the output of `get_pathway_orthologies` to only include the orthologies
#' that are related with the pathway of interest.
#'
#' @export
get_genome_orthologies <- function(data, pathway_orthologies = NULL, kegg_orthology = NULL){
  if (!"tbl_df" %in% class(data)) stop("data must be of class 'tbl_df'")
  if (!"genome_id" %in% colnames(data)) stop("data must have column named 'genome_id'")

  if (is.null(kegg_orthology)){
    kegg_orthology <- KEGGerator::kegg_orthologies
  } else {
    if (!is_kegg_tbl(kegg_orthology, "orthology")){
      stop("kegg_orthology must be a kegg_tbl with columns orthology and orthology_id", call. = FALSE)
    }
  }

  output <- data %>%
    dplyr::mutate(
      orthology_id = purrr::map(genome_id, ~{
        id <- stringr::str_replace(.x, "gn:", "")

        orths <- kegg_link_safe("orthology", id) %>%
          unique()

        if (!is.null(pathway_orthologies)) {
          orths <- orths[orths %in% pathway_orthologies]
        }

        return(orths)
      })
    ) %>%
    tidyr::unnest(orthology_id) %>%
    dplyr::left_join(kegg_orthology, by = "orthology_id")

  return(output)
}

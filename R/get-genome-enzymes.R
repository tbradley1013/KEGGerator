#' Query all enzymes related to each genome in the dataset
#'
#' @param data a tibble with genomes and genome_id or the output from `get_genome_id()`
#' @param pathway_enzymes a vector of pathway enzymes
#' @param kegg_enzyme a kegg_tbl that has the columns enzyme and enzyme_id. This
#' can be generated using the get_kegg_enzyme() function. If NULL (default) than
#' the KEGGerator::kegg_enzymes dataset will be used
#'
#' @details if `pathway_enzymes` is NULL (default) than all enzymes for each genome
#' will be returned. If a vector of enzyme_ids is provided than only the enzymes related
#' to each genome that are in this vector will be returned. This is useful if used in
#' coordination with the output of `get_pathway_enzymes` to only include the enzymes
#' that are related with the pathway of interest.
#'
#' @export
get_genome_enzymes <- function(data, pathway_enzymes = NULL, kegg_enzyme = NULL){
  if (!"tbl_df" %in% class(data)) stop("data must be of class 'tbl_df'")
  if (!"genome_id" %in% colnames(data)) stop("data must have column named 'genome_id'")

  if (is.null(kegg_enzyme)){
    kegg_enzyme <- KEGGerator::kegg_enzymes
  } else{
    if (!is_kegg_tbl(kegg_enzyme, "enzyme")){
      stop("kegg_enzyme must be a kegg_tbl with columns enzyme and enzyme_id", call. = FALSE)
    }
  }

  output <- data %>%
    dplyr::mutate(
      enzyme_id = purrr::map(genome_id, ~{
        id <- stringr::str_replace(.x, "gn:", "")

        enzymes <- kegg_link_safe("enzyme", id) %>%
          unique()

        if (!is.null(pathway_enzymes)) {
          enzymes <- enzymes[enzymes %in% pathway_enzymes]
        }

        return(enzymes)
      })
    ) %>%
    tidyr::unnest(enzyme_id) %>%
    dplyr::left_join(kegg_enzyme, by = "enzyme_id")

  return(output)
}

get_enzyme_internal <- function(orgs_id, pathway_enzymes = NULL, kegg_enzymes = NULL, verbose = FALSE, progress = TRUE){

  if (!is_orgs_id(orgs_id)) stop("orgs_id must be of class orgs_id", call. = FALSE)

  if (is.null(kegg_enzyme)){
    kegg_enzyme <- KEGGerator::kegg_enzymes
  } else{
    if (!is_kegg_tbl(kegg_enzyme, "enzyme")){
      stop("kegg_enzyme must be a kegg_tbl with columns enzyme and enzyme_id", call. = FALSE)
    }
  }

  if (!is_filtered(orgs_id)){
    warning("The orgs_id object provided has not been filtered, did you forget to run filter_orgs()?")
  }

  if (progress){
    p <- dplyr::progress_estimated(nrow(orgs_id), 10)
  }


  output <- data %>%
    dplyr::mutate(
      enzyme_id = purrr::map(genome_id, ~{
        id <- stringr::str_replace(.x, "gn:", "")

        enzymes <- kegg_link_safe("enzyme", id) %>%
          unique()

        n_hits <- length(enzymes)

        if (!is.null(pathway_enzymes)) {
          enzymes <- enzymes[enzymes %in% pathway_enzymes]
        }

        n_remain <- length(enzymes)

        if (verbose){
          cat("Linked ", crayon::red(n_hits), " (", crayon::red(n_remain),  " within the pathway specified) enzymes to genome: ", crayon::blue(.x), )
        }

        if (progress){
          p$tick()$print()
        }

        return(enzymes)
      })
    ) %>%
    tidyr::unnest(enzyme_id) %>%
    dplyr::left_join(kegg_enzyme, by = "enzyme_id")


  if (progress & verbose){
    p$stop()
  }


  return(output)

}

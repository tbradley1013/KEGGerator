#' Query all enzymes related to each genome in the dataset
#'
#' @param data a tibble with genomes and genome_id or the output from `get_genome_id()`
#' @param pathway_enzymes a vector of pathway enzymes
#' @param kegg_enzymes a kegg_tbl that has the columns enzyme and enzyme_id. This
#' can be generated using the get_kegg_enzyme() function. If NULL (default) than
#' the KEGGerator::kegg_enzymes dataset will be used
#' @param verbose logical; if TRUE the number of enzyme links for each genome
#' id will be shown as they are processed
#' @param progress logical; if TRUE (default) than a progress bar will appear if
#' the query takes longer than 10 seconds (it likely will if your data has more)
#' than only a few genomes
#'
#' @details if `pathway_enzymes` is NULL (default) than all enzymes for each genome
#' will be returned. If a vector of enzyme_ids is provided than only the enzymes related
#' to each genome that are in this vector will be returned. This is useful if used in
#' coordination with the output of `get_pathway_enzymes` to only include the enzymes
#' that are related with the pathway of interest.
#'
#' @export
get_enzymes <- function(data, pathway_enzymes, kegg_enzymes,
                        verbose, progress){
  UseMethod("get_enzymes")
}

#' @describeIn get_enzymes method for keggerator objects
#' @export
get_enzymes.keggerator <- function(data, pathway_enzymes = NULL, kegg_enzymes = NULL,
                                   verbose = FALSE, progress = TRUE){

  if (is.null(data$orgs_id)) stop("orgs_id is NULL. Have you run get_orgs_id() yet?", call. = FALSE)
  if (!is_orgs_id(data$orgs_id)) stop("object in orgs_id slot is not of class orgs_id. Did you run get_orgs_id()?")

  enzymes <- get_enzymes.orgs_id(data$orgs_id, pathway_enzymes = pathway_enzymes, kegg_enzymes = kegg_enzymes,
                                verbose = verbose, progress = progress)

  data$orgs_enzymes <- enzymes

  return(data)


}


#' @describeIn get_enzymes method for orgs_id objects
#' @export
get_enzymes.orgs_id <- function(orgs_id, pathway_enzymes = NULL, kegg_enzymes = NULL, verbose = FALSE, progress = TRUE){

  if (is.null(kegg_enzymes)){
    kegg_enzymes <- KEGGerator::kegg_enzymes
  } else{
    if (!is_kegg_tbl(kegg_enzymes, "enzyme")){
      stop("kegg_enzyme must be a kegg_tbl with columns enzyme and enzyme_id", call. = FALSE)
    }
  }

  if (!is_filtered(orgs_id)){
    warning("The orgs_id object provided has not been filtered, did you forget to run filter_orgs()?", call. = FALSE)
  }

  if (progress){
    p <- dplyr::progress_estimated(nrow(orgs_id), 10)
  }

  if (is.null(pathway_enzymes)){
    path_enzymes <- NULL
  } else if (is_keggtap(pathway_enzymes)){
    path_enzymes <- pathway_enzymes$enzyme$enzyme_id
  } else if (tibble::is_tibble(pathway_enzymes) | is.data.frame(pathway_enzymes)){
    if (!"enzyme_id" %in% colnames(pathway_enzymes)){
      stop("if pathway_enzymes is a tbl_df than it must have a column named enzyme_id")
    }

    path_enzymes <- pathway_enzymes$enzyme_id
  } else {
    path_enzymes <- pathway_enzymes
  }


  output <- orgs_id %>%
    dplyr::mutate(
      enzyme_id = purrr::map(genome_id, ~{
        id <- stringr::str_replace(.x, "gn:", "")

        enzymes <- kegg_link_safe("enzyme", id) %>%
          unique()

        n_hits <- length(enzymes)

        if (!is.null(path_enzymes)) {
          enzymes <- enzymes[enzymes %in% path_enzymes]
        }

        n_remain <- length(enzymes)

        if (verbose){
          if (!is.null(pathway_enzymes)){
            cat("Linked ", crayon::red(n_hits), " (", crayon::red(n_remain),  " within the pathway specified) enzymes to genome: ", crayon::blue(.x), "\n")
          } else {
            cat("Linked ", crayon::red(n_hits), " enzymes to genome: ", crayon::blue(.x), "\n")
          }

        }

        if (progress){
          p$tick()$print()
        }

        return(enzymes)
      })
    ) %>%
    tidyr::unnest(enzyme_id) %>%
    dplyr::left_join(kegg_enzymes, by = "enzyme_id")


  class(output) <- c("kegg_fun", class(output))
  attr(output, "query") <- "enzymes"


  if (progress & verbose){
    p$stop()
  }


  return(output)

}

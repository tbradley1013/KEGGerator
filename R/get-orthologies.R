#' Query all orthologies related to each genome in the dataset
#'
#' @param data a tibble with genomes and genome_id or the output from `get_genome_id()`
#' @param pathway_orthologies a vector of pathway orthologies
#' @param kegg_orthology a kegg_tbl with the columns orthology and orthology_id
#' which can be generated with the get_kegg_orthology function. If NULL (default)
#' the KEGGerator::kegg_orthologies dataset will be used
#' @param verbose logical; if TRUE the number of enzyme links for each genome
#' id will be shown as they are processed
#' @param progress logical; if TRUE (default) than a progress bar will appear if
#' the query takes longer than 10 seconds (it likely will if your data has more)
#' than only a few genomes
#'
#' @details if `pathway_orthologies` is NULL (default) than all orthologies for each genome
#' will be returned. If a vector of orthology_ids is provided than only the orthologies related
#' to each genome that are in this vector will be returned. This is useful if used in
#' coordination with the output of `get_pathway_orthologies` to only include the orthologies
#' that are related with the pathway of interest.
#'
#' @export
get_orthologies <- function(data, pathway_orthologies, kegg_orthology, verbose, progress){
  UseMethod("get_orthologies")
}

#' @describeIn get_orthologies method for keggerator objects
#' @export
get_orthologies.keggerator <- function(data, pathway_orthologies = NULL, kegg_orthology = NULL,
                                       verbose = FALSE, progress = TRUE){
  if (is.null(data$orgs_id)) stop("orgs_id is NULL. Have you run get_orgs_id() yet?", call. = FALSE)
  if (!is_orgs_id(data$orgs_id)) stop("object in orgs_id slot is not of class orgs_id. Did you run get_orgs_id()?")

  orths <- get_orthologies.orgs_id(data$orgs_id, pathway_orthologies = pathway_orthologies,
                                   kegg_orthology = kegg_orthology, verbose = verbose, progress = progress)

  data$orgs_orthologies <- orths

  return(data)
}


#' @describeIn get_orthologies method for orgs_id objects
#' @export
get_orthologies.orgs_id <- function(orgs_id, pathway_orthologies, kegg_orthology = NULL,
                                    verbose = FALSE, progress = TRUE){

  if (is.null(kegg_orthology)){
    kegg_orthology <- KEGGerator::kegg_orthologies
  } else{
    if (!is_kegg_tbl(kegg_orthology, "enzyme")){
      stop("kegg_enzyme must be a kegg_tbl with columns enzyme and enzyme_id", call. = FALSE)
    }
  }

  if (!is_filtered(orgs_id)){
    warning("The orgs_id object provided has not been filtered, did you forget to run filter_orgs()?", call. = FALSE)
  }

  if (progress){
    p <- dplyr::progress_estimated(nrow(orgs_id), 10)
  }

  output <- data %>%
    dplyr::mutate(
      orthology_id = purrr::map(genome_id, ~{
        id <- stringr::str_replace(.x, "gn:", "")

        orths <- kegg_link_safe("orthology", id) %>%
          unique()

        n_hits <- length(orths)

        if (!is.null(pathway_orthologies)) {
          orths <- orths[orths %in% pathway_orthologies]
        }

        n_remain <- length(orths)

        if (verbose){
          if (!is.null(pathway_orthologies)){
            cat("Linked ", crayon::red(n_hits), " (", crayon::red(n_remain),  " within the pathway specified) orthologies to genome: ", crayon::blue(.x), "\n")
          } else {
            cat("Linked ", crayon::red(n_hits), " orthologies to genome: ", crayon::blue(.x), "\n")
          }
        }

        if (progress){
          p$tick()$print()
        }

        return(orths)
      })
    ) %>%
    tidyr::unnest(orthology_id) %>%
    dplyr::left_join(kegg_orthology, by = "orthology_id")

  class(output) <- c("kegg_fun", class(output))
  attr(output, "query") <- "orthologies"


  if (progress & verbose){
    p$stop()
  }


  return(output)

}


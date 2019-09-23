#' retrieve KEGG genome id for all organims
#'
#' @param data dataset of genome names to get id for
#'
#' @details This took ~4.4 minutes to run for 400 organisms
#'
#' @export
get_org_ids <- function(data, verbose, progress){
  UseMethod("get_org_ids")
}

get_org_ids.orgs_tbl <- function(data, verbose = TRUE, progress = TRUE){
  if (!is_orgs_tbl(data)) stop("data must have column named genome", call. = FALSE)

  unq_orgs <- dplyr::distinct(data, genome)

  if (progress){
    p <- dplyr::progress_estimated(nrow(unq_orgs), 10)
  }

  org_hits <- unq_orgs %>%
    dplyr::mutate(genome_query = purrr::map(genome, ~{
      query <- kegg_find_safe("genome", .x)

      if (verbose){
        cat(
          "Finding genomes that match ", crayon::red(.x), ": ",
          crayon::blue(length(query)), "\n", sep = ""
        )
      }

      if (progress){
        p$tick()$print()
      }

      return(query)
    })) %>%
    dplyr::filter(purrr::map(genome_query, length) > 0,
                  !is.na(genome_query)) %>%
    dplyr::mutate(
      genome_id = purrr::map(genome_query, ~names(.x))
    ) %>%
    tidyr::unnest(genome_query, genome_id) %>%
    dplyr::distinct() %>%
    tidyr::separate(genome_query, into = c("genome_name", "genome_desc"), sep = "; ")


  kegg_uncert <- org_hits %>%
    dplyr::group_by(genome) %>%
    dplyr::summarise(
      n_hits = sum(!is.na(genome_id))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      uncert = dplyr::if_else(n_hits == 0, 1, 1-(1/n_hits))
    ) %>%
    dplyr::inner_join(ps_kegg$orgs_tbl, by = "genome") %>%
    dplyr::arrange(otu_id)

  class(org_hits) <- c("orgs_id", class(org_hits))
  attr(orgs_hits, "filtered") <- FALSE
  class(kegg_uncert) <- c("uncert_tbl", class(kegg_uncert))

  output <- list(orgs_id = org_hits, kegg_uncert = kegg_uncert)
  class(output) <- c("orgs_list", class(output))

  if (verbose & progress){
    p$stop()
  }

  return(output)

}

#' @describeIn get_org_ids method for keggerator
#' @export
get_org_ids.keggerator <- function(data, verbose = FALSE, progress = TRUE){

  orgs <- data$orgs_tbl

  ids <- get_org_ids.orgs_tbl(orgs, verbose = verbose, progress = progress)

  data$orgs_id <- ids$orgs_id
  data$kegg_uncert <- ids$kegg_uncert

  # calculating the total uncertainty
  if (!is.null(data$species_uncert) & !is.null(data$kegg_uncert)){
    data <- keggerator_uncertainty(data)
  }

  return(data)

}

#' @describeIn get_org_ids method for orgs_list
#' @export
get_org_ids.orgs_list <- function(data, verbose = FALSE, progress = TRUE){

  orgs <- data$orgs_tbl

  ids <- get_org_ids.orgs_tbl(orgs, verbose = verbose, progress = progress)

  data$orgs_id <- ids$orgs_id
  data$kegg_uncert <- ids$kegg_uncert

  # calculating the total uncertainty
  if (!is.null(data$species_uncert) & !is.null(data$kegg_uncert)){
    data <- keggerator_uncertainty(data)
  }

  return(data)

}


#' @export
is_orgs_id <- function(x){
  all(c(inherits(x, "orgs_id"), "genome_id" %in% colnames(x)))
}

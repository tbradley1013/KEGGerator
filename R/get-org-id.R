#' retrieve KEGG genome id for all organims
#'
#' @param data dataset of genome names to get id for
#'
#' @details This took ~4.4 minutes to run for 400 organisms
#'
#' @export
get_org_ids <- function(data, verbose){
  UseMethod("get_org_ids")
}

get_org_ids.orgs_tbl <- function(data, verbose = TRUE){
  if (!is_orgs_tbl) stop("data must have column named genome", call. = FALSE)

  unq_orgs <- dplyr::distinct(data, genome)

  output <- unq_orgs %>%
    dplyr::mutate(genome_query = purrr::map(genome, ~{
      query <- kegg_find_safe("genome", .x)

      if (verbose){
        cat(
          "Matching Genomes for ", crayon::red(.x), ": ",
          crayon::blue(length(query)), "\n", sep = ""
        )
      }

      return(query)
    })) %>%
    dplyr::filter(purrr::map(genome_query, length) > 0,
                  !is.na(genome_query)) %>%
    dplyr::mutate(
      genome_id = purrr::map(genome_query, ~names(.x))
      # n_genomes = purrr::map(genome_query, length)
    ) %>%
    tidyr::unnest(genome_query, genome_id) %>%
    dplyr::distinct()





  return(output)

}

#' @export
get_genome_id <- function(data, verbose = FALSE) {
  if (any(class(data) == "phyloseq")) data <- genomes_tibble(data)

  if (!"tbl_df" %in% class(data)) stop("data must be either of class tbl_df or phyloseq")
  if (!"genome" %in% colnames(data)) stop("there must be a column named 'genome' in data")

  output <- data %>%
    dplyr::mutate(genome_query = purrr::map(genome, ~{
      query <- kegg_find_safe("genome", .x)

      if (verbose){
        cat(
          "Matching Genomes for ", crayon::red(.x), ": ",
          crayon::blue(length(query)), "\n", sep = ""
        )
      }

      return(query)
    })) %>%
    dplyr::filter(purrr::map(genome_query, length) > 0,
                  !is.na(genome_query)) %>%
    dplyr::mutate(genome_id = purrr::map(genome_query, ~names(.x))) %>%
    tidyr::unnest(genome_query, genome_id) %>%
    dplyr::distinct()





  return(output)
}

#' retrieve KEGG genome id for all genomes
#'
#' @param data dataset of genome names to get id for
#'
#' @importFrom magrittr %>%
#'
#' @details This took ~4.4 minutes to run for 400 genomes
#'
#' @export
get_genome_id <- function(data) {
  if (class(data) == "phyloseq") data <- tibble_genomes(data)

  if (!"tbl_df" %in% class(data)) stop("data must be either of class tbl_df or phyloseq")
  if (!"genome" %in% colnames(data)) stop("there must be a column named 'genome' in data")

  output <- data %>%
    dplyr::mutate(genome_query = purrr::map(genome, ~KEGGREST::keggFind("genome", .x))) %>%
    dplyr::filter(purrr::map(genome_query, length) > 0) %>%
    dplyr::mutate(genome_id = purrr::map(genome_query, ~names(.x))) %>%
    tidyr::unnest(genome_query, genome_id) %>%
    dplyr::distinct()

  return(output)
}

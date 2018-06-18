#' extract genomes from tax table or tax tibble
#'
#' @param data either a phyloseq object or the output of tibble_tax
#'
#' @importFrom magrittr %>%
#'
#' @export
tibble_genomes <- function(data){
  if (class(data) == "phyloseq") data <- tibble_tax(data)

  if (!"tbl_df" %in% class(data)) stop("data must be either of class tbl_df or phyloseq")

  taxa <- data %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::filter(!is.na(Genus)) %>%
    dplyr::mutate(Species = dplyr::if_else(is.na(Species), "", Species))

  max_uncertainty <- taxa$Species %>%
    stringr::str_count("\\/") %>%
    max()

  taxa_cols <- purrr::map_chr(1:(max_uncertainty + 1), ~paste0("genome", .x))

  output <- taxa %>%
  {suppressWarnings(tidyr::separate(., Species, into = taxa_cols))} %>%
    tidyr::gather(key = key, value = Species, !!dplyr::quo(taxa_cols)) %>%
    dplyr::filter(!is.na(Species)) %>%
    tidyr::unite("genome", Genus, Species, sep = " ") %>%
    dplyr::mutate(genome = stringr::str_trim(genome)) %>%
    dplyr::select(genome)

  return(output)
}

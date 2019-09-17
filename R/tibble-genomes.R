#' extract genomes from tax table or tax tibble
#'
#' @param data either a phyloseq object or the output of tibble_tax
#' @param drop_taxa logical; should the taxonomy be removed from the output
#' @param strict logical; if set to TRUE, it will not keep any taxa that have
#' not been classified to the species level.
#'
#' @export
genomes_tibble <- function(data, drop_taxa, strict = FALSE){
  UseMethod("genomes_tibble")
}



#' @export
genomes_tibble <- function(data, drop_taxa = TRUE, strict = FALSE){
  if (any(class(data) == "phyloseq")) data <- tax_tibble(data)

  if (!"tbl_df" %in% class(data)) stop("data must be either of class tbl_df or phyloseq")

  taxa <- data %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::filter(!is.na(Genus)) %>%
    dplyr::mutate(Species = dplyr::if_else(is.na(Species), "", Species))

  # removing the taxa that do not have any species assigned
  if (strict){
    taxa <- dplyr::filter(taxa, Species != "")
  }



  max_uncertainty <- taxa$Species %>%
    stringr::str_count("\\/") %>%
    max()

  taxa_cols <- purrr::map_chr(1:(max_uncertainty + 1), ~paste0("genome", .x))

  output <- taxa %>%
    {suppressWarnings(tidyr::separate(., Species, into = taxa_cols))} %>%
    tidyr::gather(key = key, value = Species, !!dplyr::quo(taxa_cols)) %>%
    dplyr::filter(!is.na(Species)) %>%
    tidyr::unite("genome", Genus, Species, sep = " ") %>%
    dplyr::mutate(genome = stringr::str_trim(genome))

  if (drop_taxa) output <- dplyr::select(output, otu, genome)

  return(output)
}


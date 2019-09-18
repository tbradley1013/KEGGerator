#' extract organisms from tax table or tax tibble
#'
#' @param data either a phyloseq object or the output of tibble_tax
#' @param drop_taxa logical; should the taxonomy be removed from the output
#' @param strict logical; if set to TRUE, it will not keep any taxa that have
#' not been classified to the species level.
#' @param sep what is used to separate multiple species in the species column.
#' Must be valid regex targeting the seperator. Default is "\\/" to match the
#' default separator used during taxonomic classification by \code{\link[dada2]{addSpecies}}
#'
#' @details
#' When an object of class tax_tbl is passed to orgs_tibble, a list with two
#' tibbles will be returned. The first tibble is the orgs tibble and the second
#' is the uncertainty tibble which is comprised of the percent uncertainty [0-1]
#' of the species level assignment of each otu. If the otu was not assigned to
#' the species level than the uncertainty level is 1. If the otu was assigned to
#' only a single species than the uncertainty level is 0. If the species level
#' was assigned to N possible species, than the uncertainty level is 1-(1/N).
#'
#' If an object of class keggerator is passed to orgs_tibble than the same two
#' tbls are returned, but rather than being in a list alone, they are added to
#' the keggerator object that is given in the orgs_tbl and species_uncert
#' slots, respectively.
#'
#' @export
orgs_tibble <- function(data, drop_taxa, strict, sep){
  UseMethod("orgs_tibble")
}

#' @describeIn orgs_tibble method for tax_tbl
#' @export
orgs_tibble.tax_tbl <- function(data, drop_taxa = TRUE, strict = FALSE, sep = "\\/"){
  if (!ids_match(data)){
    warning("tax_tbl provided has not had otu_ids verified against otu_tbl and may cause issues later if they do not match.",
            "Run KEGGerator:::check_otu_id() on your phyloseq data to ensure OTUs match.", call. = FALSE)
  }

  # removing all otus that are not defined to at least the genus level
  taxa <- data %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::filter(!is.na(Genus)) %>%
    dplyr::mutate(Species = dplyr::if_else(is.na(Species), "", Species))

  # getting all otu_ids that are not defined at the species level
  spec_undef <- taxa$otu_id[taxa$Species == ""]

  # removing the taxa that do not have any species assigned
  if (strict){
    taxa <- dplyr::filter(taxa, Species != "")
  }

  # finding the max uncertainty in any of the species columns
  max_uncertainty <- taxa$Species %>%
    stringr::str_count(sep) %>%
    max()

  taxa_cols <- purrr::map_chr(1:(max_uncertainty + 1), ~paste0("genome", .x))

  orgs <- taxa %>%
    {suppressWarnings(tidyr::separate(., Species, into = taxa_cols, sep = sep))} %>%
    tidyr::gather(key = key, value = Species, !!dplyr::quo(taxa_cols)) %>%
    dplyr::filter(!is.na(Species)) %>%
    tidyr::unite("genome", Genus, Species, sep = " ") %>%
    dplyr::mutate(genome = stringr::str_trim(genome))

  if (drop_taxa) orgs <- dplyr::select(orgs, otu_id, genome)

  otu_species_uncert <- orgs %>%
    dplyr::group_by(otu_id) %>%
    dplyr::summarize(n_spec = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(uncert = dplyr::if_else(otu_id %in% spec_undef, 1, 1-(1/n_spec)))

  class(orgs) <- c("orgs_tbl", class(orgs))
  class(otu_species_uncert) <- c("uncert_tbl", class(otu_species_uncert))

  output <- list(orgs_tibble = orgs, species_uncert = otu_species_uncert)

  class(output) <- c("orgs_list", class(output))

  return(output)
}

#' @describeIn orgs_tibble method for keggerator
#' @export
orgs_tibble.keggerator <- function(data, drop_taxa = TRUE, strict = FALSE, sep = "\\/"){

  tax <- data$tax_tbl

  orgs <- orgs_tibble.tax_tbl(tax, drop_taxa = drop_taxa, strict = strict, sep = sep)

  data$orgs_tbl <- orgs$orgs_tibble
  data$species_uncert <- orgs$species_uncert

  return(data)

}

#' @describeIn orgs_tibble method for phyloseq object
#' @export
orgs_tibble.phyloseq <- function(data, drop_taxa = TRUE, strict = FALSE, sep = "\\/"){

  data <- as_keggerator(data)

  output <- orgs_tibble.keggerator(data, drop_taxa = drop_taxa, strict = strict, sep = sep)

  return(output)
}

is_orgs_list <- function(x){
  inherits(x, "orgs_list")
}

is_orgs_tbl <- function(x){
  all(c(inherits(x, "orgs_tbl"), "genome" %in% colnames(x)))
}

is_uncert_tbl <- function(x){
  inherits(x, "uncert_tbl")
}

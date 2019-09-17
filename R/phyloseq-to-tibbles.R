# These functions convert the phyloseq into separate tibbles based on which table
# you are interested in. These will likely be used inside other functions to
# get the data in the format I like. But they could be useful so could be
# exported when we get that far
#===============================================================================

#' Convert OTU table to tibble
#'
#' @param data a phyloseq object
#' @export
otu_tibble <- function(data){
  UseMethod("otu_tibble")
}


#' @export
otu_tibble.phyloseq <- function(data) {
  # if (class(data) != "phyloseq") stop("data must be a phyloseq object")

  otu <- phyloseq::otu_table(data)
  if (!phyloseq::taxa_are_rows(otu)) {
    otu <- t(otu)
  }

  output <- otu %>%
    as.data.frame() %>%
    tibble::rownames_to_column("otu") %>%
    tibble::as_tibble()

  return(output)
}

#' Convert tax table to tibble
#'
#' @param data a phyloseq object
#' @export
tax_tibble <- function(data){
  UseMethod("tax_tibble")
}


#' @export
tax_tibble.phyloseq <- function(data){
  if (class(data) != "phyloseq") stop("data must be a phyloseq object")

  tax <- phyloseq::tax_table(data)

  output <- tax %>%
    as.data.frame() %>%
    tibble::rownames_to_column("otu") %>%
    tibble::as_tibble()

  return(output)
}


#' Convert sam_data to tibble
#'
#' @param data a phyloseq object
#'
#' @export
sam_tibble <- function(data){
  UseMethod("sam_tibble")
}


#' @export
sam_tibble.phyloseq <- function(data){
  sam_data <- phyloseq::sam_data(data)

  output <- sam_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble()

  return(output)
}

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

  otu <- phyloseq::otu_table(data)
  if (!phyloseq::taxa_are_rows(otu)) {
    otu <- t(otu)
  }
  id_match <- check_otu_id(data)

  output <- otu %>%
    as.data.frame() %>%
    tibble::rownames_to_column("otu") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(otu_id = dplyr::row_number()) %>%
    dplyr::select(otu_id, otu, dplyr::everything())

  class(output) <- c("tax_tbl", class(output))
  attr(output, "id_match") <- id_match

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

  tax <- phyloseq::tax_table(data)
  id_match <- check_otu_id(data)

  output <- tax %>%
    as.data.frame() %>%
    tibble::rownames_to_column("otu") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(otu_id = dplyr::row_number()) %>%
    dplyr::select(otu_id, otu, dplyr::everything())

  class(output) <- c("otu_tbl", class(output))
  attr(output, "id_match") <- id_match

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

  class(output) <- c("sam_tbl", class(output))

  return(output)
}


#' Create a otu reference table
#'
#' @param data an object to convert to otu_ref
#'
#' @export
otu_ref <- function(data){
  UseMethod("otu_ref")
}


#' @export
otu_ref.phyloseq <- function(data){
  ids_match <- check_otu_id(data)
  if (!ids_match) stop("The OTU names do not match between the tax_table and otu_table in phyloseq object provided", call. = FALSE)

  otus <- rownames(tax_table(data))

  output <- tibble::tibble(
    otu_id = seq_along(otus),
    otu = otus
  )

  class(output) <- c("otu_ref", class(output))
  attr(output, "id_match") <- ids_match

  return(output)
}


check_otu_id <- function(data){
  tax <- phyloseq::tax_table(data)
  otu <- phyloseq::otu_table(data)
  if (!phyloseq::taxa_are_rows(otu)) {
    otu <- t(otu)
  }

  identical(rownames(tax), rownames(otu))
}

ids_match <- function(x){
  attr(x, "id_match")
}


is_keggerator <- function(x){
  inherits(x, "keggerator")
}

is_otu_tbl <- function(x){
  inherits(x, "otu_tbl")
}

is_tax_tbl <- function(x){
  inherits(x, "tax_tbl")
}

is_sam_tbl <- function(x){
  inherits(x, "sam_tbl")
}

is_otu_ref <- function(x){
  inherits(x, "otu_ref")
}

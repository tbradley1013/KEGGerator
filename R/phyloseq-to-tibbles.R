# These functions convert the phyloseq into separate tibbles based on which table
# you are interested in. These will likely be used inside other functions to
# get the data in the format I like. But they could be useful so could be
# exported when we get that far
#===============================================================================

#' Convert OTU table to tibble
#'
#' @param data a phyloseq object
#' @param quiet whether to pass warnings about otu matching or not
#' @export
otu_tibble <- function(data, quiet){
  UseMethod("otu_tibble")
}


#' @describeIn otu_tibble
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

  class(output) <- c("otu_tbl", class(output))
  attr(output, "id_match") <- id_match

  return(output)
}

#' @describeIn otu_tibble
#' @export
otu_tibble.otu_table <- function(data, quiet = FALSE){
  otu <- data

  if (!quiet) warning("Passing an otu_table directly to otu_tibble does not allow for verification that OTU names and order match between otu_table and tax_table", call. = FALSE)

  if (!phyloseq::taxa_are_rows(otu)) {
    otu <- t(otu)
  }


  output <- otu %>%
    as.data.frame() %>%
    tibble::rownames_to_column("otu") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(otu_id = dplyr::row_number()) %>%
    dplyr::select(otu_id, otu, dplyr::everything())

  class(output) <- c("otu_tbl", class(output))
  attr(output, "id_match") <- NULL

  return(output)
}

#' Convert tax table to tibble
#'
#' @param data a phyloseq object
#' @param quiet whether to pass warnings about otu matching or not
#' @export
tax_tibble <- function(data, quiet){
  UseMethod("tax_tibble")
}


#' @describeIn tax_tibble
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

  class(output) <- c("tax_tbl", class(output))
  attr(output, "id_match") <- id_match

  return(output)
}

#' @describeIn tax_tibble
#' @export
tax_tibble.taxonomyTable <- function(data, quiet = FALSE){
  tax <- data

  if (!quiet) warning("Passing an otu_table directly to otu_tibble does not allow for verification that OTU names and order match between otu_table and tax_table", call. = FALSE)


  output <- tax %>%
    as.data.frame() %>%
    tibble::rownames_to_column("otu") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(otu_id = dplyr::row_number()) %>%
    dplyr::select(otu_id, otu, dplyr::everything())

  class(output) <- c("tax_tbl", class(output))
  attr(output, "id_match") <- NULL

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


#' @describeIn sam_tibble
#' @export
sam_tibble.phyloseq <- function(data){
  sam_data <- phyloseq::sample_data(data)

  output <- sam_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble() %>%
    dplyr::select(sample, dplyr::everything())

  class(output) <- c("sam_tbl", class(output))

  return(output)
}

#' @describeIn sam_tibble
#' @export
sam_tibble.sample_data <- function(data){

  output <- sam_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble()  %>%
    dplyr::select(sample, dplyr::everything())

  class(output) <- c("sam_tbl", class(output))

  return(output)

}


#' Create a otu reference table
#'
#' @param data an object to convert to otu_ref
#'
#' @export
otu_ref <- function(data, otu, tax){
  UseMethod("otu_ref")
}


#' @describeIn otu_ref
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

#' @describeIn otu_ref
#' @export
otu_ref.otu_table <- function(otu, tax){
  ids_match <- check_otu_id(otu, tax)
  if (!ids_match) stop("The OTU names do not match between the tax_table and otu_table in the otu_table and taxonomyTable objects provided", call. = FALSE)

  otus <- rownames(tax)

  output <- tibble::tibble(
    otu_id = seq_along(otus),
    otu = otus
  )

  class(output) <- c("otu_ref", class(output))
  attr(output, "id_match") <- ids_match

  return(output)
}

#' @describeIn otu_ref
#' @export
otu_ref.taxonomyTable <- otu_ref.otu_table

#' @describeIn otu_ref
#' @export
otu_ref.tax_tbl <- function(tax){
  id_match <- attr(tax, "id_match")
  if (is.null(id_match) | !id_match) stop("tax_tbl provided has not been checked for otu agreement between tax table and otu table", call. = FALSE)

  output <- dplyr::select(tax, otu_id, otu)

  prev_class <- class(output)
  class(output) <- c("otu_ref", prev_class[prev_class != "tax_tbl"])
  return(output)
}

#' @describeIn otu_ref
#' @export
otu_ref.otu_tbl <- function(otu){
  id_match <- attr(otu, "id_match")
  if (is.null(id_match) | !id_match) stop("otu_tbl provided has not been checked for otu agreement between tax table and otu table", call. = FALSE)

  output <- dplyr::select(otu, otu_id, otu)

  prev_class <- class(output)
  class(output) <- c("otu_ref", prev_class[prev_class != "otu_tbl"])
  return(output)
}

# This function checks whether the otu names match for
# both the taxonomy table and the otu table either within
# a single phyloseq object or in separately provided taxonomyTable
# and otu_table objects
check_otu_id <- function(data, otu, tax){
  UseMethod("check_otu_id")
}

check_otu_id.phyloseq <- function(data){
  tax <- phyloseq::tax_table(data)
  otu <- phyloseq::otu_table(data)
  if (!phyloseq::taxa_are_rows(otu)) {
    otu <- t(otu)
  }

  identical(rownames(tax), rownames(otu))
}


check_otu_id.otu_table <- function(otu, tax){
  if (!inherits(otu, "otu_table")) stop("otu must be of class otu_table")
  if (!inherits(tax, "taxonomyTable")) stop("tax must be of class tax_table")
  if (!phyloseq::taxa_are_rows(otu)){
    otu <- t(otu)
  }

  identical(rownames(tax), rownames(otu))
}

check_otu_id.taxonomyTable <- check_otu_id.otu_table


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

#' A keggerator object
#'
#' A class of object that contains all of the
#' required components for analysis with
#' KEGGerator
#'
#' @return an object of class keggerator
#'
#' @export
keggerator <- function(tax_tbl = NULL, otu_tbl = NULL, sam_tbl = NULL,
                       otu_ref = NULL, orgs_tbl = NULL, species_uncert = NULL, orgs_id = NULL,
                       orgs_filt = NULL, orgs_enzymes = NULL, orgs_orthologies = NULL,
                       orgs_genes = NULL){
  # These first four arguments are required in order to build the keggerator object
  if (!null_check_req(tax_tbl, is_tax_tbl)) stop("argument tax_tbl must not be null and be of class tax_tbl", call. = FALSE)
  if (!null_check_req(otu_tbl, is_otu_tbl)) stop("argument otu_tbl must not be null and be of class otu_tbl", call. = FALSE)
  if (!null_check_req(sam_tbl, is_sam_tbl)) stop("argument sam_tbl must not be null and be of class sam_tbl", call. = FALSE)
  if (!null_check_req(otu_ref, is_otu_ref)) stop("argument otu_ref must not be null and be of class otu_ref", call. = FALSE)

  # Checking if the id_match attribute is set for the tax_tbl and otu_tbl
  # arguments. If it is not set, than it will be set using the id_match
  # attribute from the otu_ref object. This is done because if the user
  # passes otu_table and taxonomyTable objects to create the otu_tbl and
  # tax_tbl objects, respectively, this attribute can not be set.
  if (is.null(attr(tax_tbl, "id_match"))){
    attr(tax_tbl, "id_match") <- attr(otu_ref, "id_match")
  }

  if (is.null(attr(otu_tbl, "id_match"))){
    attr(otu_tbl, "id_match") <- attr(otu_ref, "id_match")
  }

  output <- structure(
    list(
      "tax_tbl" = tax_tbl,
      "otu_tbl" = otu_tbl,
      "sam_tbl" = sam_tbl,
      "otu_ref" = otu_ref,
      "orgs_tbl" = orgs_tbl,
      "species_uncert" = species_uncert,
      "orgs_id" = orgs_id,
      "orgs_filt" = orgs_filt,
      "orgs_enzymes" = orgs_enzymes,
      "orgs_orthologies" = orgs_orthologies,
      "orgs_genes" = orgs_genes
    ),
    class = "keggerator"
  )

  return(output)
}

null_check_req <- function(x, .f){
  if (is.null(x)) return(FALSE)

  .f(x)
}

null_check_opt <- function(x, .f) {
  if (is.null(x)) return(FALSE)

  .f(x)
}

#' Convert an object into class keggerator
#'
#' @param ps a phyloseq object to be converted to a keggerator object
#' @param ... additional arguements to be passed to \code{\link[KEGGerator]{keggerator}}
#'
#' @export
as_keggerator <- function(ps, ...){
  UseMethod("as_keggerator")
}

#' @describeIn as_keggerator method for phyloseq objects
#' @export
as_keggerator.phyloseq <- function(ps, ...){
  # browser()

  otu <- otu_tibble(ps)
  tax <- tax_tibble(ps)
  sam <- sam_tibble(ps)
  otu_ref <- otu_ref(ps)

  output <- keggerator(tax_tbl = tax, otu_tbl = otu, sam_tbl = sam, otu_ref = otu_ref, ...)

  return(output)
}

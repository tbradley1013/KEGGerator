#' A keggerator object
#'
#' A class of object that contains all of the
#' required components for analysis with
#' KEGGerator
#'
keggerator <- function(tax_tbl = NULL, otu_tbl = NULL, sam_tbl = NULL,
                       otu_ref = NULL, orgs_tbl = NULL, orgs_id = NULL,
                       orgs_filt = NULL, orgs_enzymes = NULL, orgs_orthologies = NULL,
                       orgs_genes = NULL){
  if (!null_check(tax_tbl, is_tax_tbl)) stop("argument tax_tbl must be of class tax_tbl")
  if (!null_check(otu_tbl, is_otu_tbl)) stop("argument otu_tbl must be of class otu_tbl")
  if (!null_check(sam_tbl, is_sam_tbl)) stop("argument sam_tbl must be of class sam_tbl")

  output <- structure(
    list(
      "tax_tbl" = tax_tbl,
      "otu_tbl" = otu_tbl,
      "sam_tbl" = sam_tbl,
      "otu_ref" = otu_ref,
      "orgs_tbl" = orgs_tbl,
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

null_check <- function(x, .f){
  if (is.null(x)) return(TRUE)

  .f(x)
}

as_keggerator <- function(x){
  UseMethod("as_keggerator")
}

as_keggerator.phyloseq <- function(x){

}

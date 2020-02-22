#' A keggerator object
#'
#' A class of object that contains all of the
#' required components for analysis with
#' KEGGerator
#'
#' @return an object of class keggerator
#'
#' @export
keggerator <- R6::R6Class(
  "keggerator",
  public = list(
    ps = NULL,
    tax_tbl = NULL,
    otu_tbl = NULL,
    sam_tbl = NULL,
    otu_ref = NULL,
    orgs_tbl = NULL,
    species_uncert = NULL,
    orgs_id = NULL,
    kegg_uncert = NULL,
    total_uncert = NULL,
    orgs_filt = NULL,
    orgs_enzymes = NULL,
    orgs_orthologies = NULL,

    initialize = function(ps){
      stopifnot(inherits(ps, "phyloseq"))
      cat("Intializing KEGGerator...\n")
      tictoc::tic("Finished Converting phyloseq object")
      self$ps <- ps
      cat("Converting Tax Table\n")
      tictoc::tic("Finished Converting Tax Table")
      self$tax_tbl <- tax_tibble(ps)
      tictoc::toc()
      cat("Converting OTU table\n")
      tictoc::tic("Finished Converting OTU Table")
      self$otu_tbl <- otu_tibble(ps)
      tictoc::toc()
      cat("Converting Sample Table\n")
      tictoc::tic("Finished Converting Sample Table")
      self$sam_tbl <- sam_tibble(ps)
      tictoc::toc()
      cat("Creating OTU Reference\n")
      tictoc::tic("Finished Creating OTU Reference")
      self$otu_ref <- otu_ref(ps)
      tictoc::toc()
      tictoc::toc()
    },


  )
)
null_check_req <- function(x, .f){
  if (is.null(x)) return(FALSE)

  .f(x)
}

null_check_opt <- function(x, .f) {
  if (is.null(x)) return(TRUE)

  .f(x)
}

verbose <- function(msg, expr){
  cat(glue::glue("{msg}"), "\n")
  tictoc::tic(glue::glue("Finished {msg}"))
  fn <- exprToFunction(expr)
  out <- fn()
  tictoc::toc()
  return(out)
}

# From Shiny
makeFunction <- function(args = pairlist(), body, env = parent.frame()) {
  eval(call("function", args, body), env)
}

exprToFunction <- function(expr, env=parent.frame(), quoted=FALSE) {
  if (!quoted) {
    expr <- eval(substitute(substitute(expr)), parent.frame())
  }

  # expr is a quoted expression
  makeFunction(body=expr, env=env)
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
  output <- keggerator@new(ps)
  return(output)
}




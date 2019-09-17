#===============================================================================
# Functions to download the latest KEGG datasets
#
# Tyler Bradley
# 2019-09-14
#===============================================================================

#' Get Updated KEGG datasets
#'
#' @name get_datasets
#'
#' @details
#' These functions will allow users to get updated datasets rather if they do
#' not want to use the ones supplied in this package. While I will try to keep
#' these datasets updated, it will allow for the flexibility for users to
#' have the most updated data in KEGG
NULL

#' @rdname get_datasets
#' @export
get_kegg_pathway <- function(){
  kegg_pathways_list <- kegg_list_safe("pathway")

  if (is.na(kegg_pathways_list)) stop("Download of KEGG pathway information failed", call. = FALSE)

  pathways <- tibble::tibble(
    pathway_id = names(kegg_pathways_list),
    pathway = kegg_pathways_list
  )

  class(pathways) <- c("kegg_tbl", class(pathways))

  return(pathways)
}


#' @rdname get_datasets
#' @export
get_kegg_module <- function(){
  kegg_module_list <- kegg_list_safe("module")

  if (is.na(kegg_module_list)) stop("Download of KEGG module information failed", call. = FALSE)

  modules <- tibble::tibble(
    module_id = names(kegg_module_list),
    module = kegg_module_list
  )

  class(modules) <- c("kegg_tbl", class(modules))

  return(modules)
}

#' @rdname get_datasets
#' @export
get_kegg_enzyme <- function(){
  kegg_enzyme_list <- kegg_list_safe("enzyme")

  if (is.na(kegg_enzyme_list)) stop("Download of KEGG enzyme information failed", call. = FALSE)

  enzymes <- tibble::tibble(
    enzyme_id = names(kegg_enzyme_list),
    enzyme = kegg_enzyme_list
  )

  class(enzymes) <- c("kegg_tbl", class(enzymes))

  return(enzymes)
}

#' @rdname get_datasets
#' @export
get_kegg_orthology <- function(){
  kegg_orthology_list <- kegg_list_safe("orthology")

  if (is.na(kegg_orthology_list)) stop("Download of KEGG orthology information failed", call. = FALSE)

  orthologies <- tibble::tibble(
    orthology_id = names(kegg_orthology_list),
    orthology = kegg_orthology_list
  )

  class(orthologies) <- c("kegg_tbl", class(orthologies))

  return(orthologies)
}

#' Check whether KEGGerator datasets have the correct class
#'
#' @export
is_kegg_tbl <- function(x, type){
  cols <- paste0(type, c("_id", ""))

  if (inherits(x, "kegg_tbl")){
    if (all(cols %in% colnames(x))){
      return(TRUE)
    }
  }

  return(FALSE)
}

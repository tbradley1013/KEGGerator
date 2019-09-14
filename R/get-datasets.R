#===============================================================================
# Functions to download the latest KEGG datasets
#
# Tyler Bradley
# 2019-09-14
#===============================================================================

#' Get Updated KEGG datasets
#' @name get_datasets
NULL

#' @rdname get_datasets
get_kegg_pathway <- function(){
  kegg_pathways_list <- kegg_list_safe("pathway")

  if (is.na(kegg_pathways_list)) stop("Download of KEGG pathway information failed", call. = FALSE)

  kegg_pathways <- tibble::tibble(
    pathway_id = names(kegg_pathways_list),
    pathway = kegg_pathways_list
  )

  return(kegg_pathways)
}

#' @rdname get_datasets
get_kegg_module <- function(){
  kegg_module_list <- kegg_list_safe("module")

  if (is.na(kegg_module_list)) stop("Download of KEGG module information failed", call. = FALSE)

  kegg_modules <- tibble::tibble(
    module_id = names(kegg_module_list),
    module = kegg_module_list
  )
}

#' @rdname get_datasets
get_kegg_enzyme <- function(){
  kegg_enzyme_list <- kegg_list_safe("enzyme")

  if (is.na(kegg_enzyme_list)) stop("Download of KEGG enzyme information failed", call. = FALSE)

  kegg_enzymes <- tibble::tibble(
    enzyme_id = names(kegg_enzyme_list),
    enzyme = kegg_enzyme_list
  )
}

#' @rdname get_datasets
get_kegg_orthology <- function(){
  kegg_orthology_list <- kegg_list_safe("orthology")

  if (is.na(kegg_enzyme_list)) stop("Download of KEGG orthology information failed", call. = FALSE)

  kegg_orthologies <- tibble::tibble(
    orthology_id = names(kegg_orthology_list),
    orthology = kegg_orthology_list
  )
}

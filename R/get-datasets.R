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

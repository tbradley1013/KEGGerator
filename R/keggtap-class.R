# function to create keggtap object

#' Create a keggtap object
#'
#' @param pathway_name a character or character vector that specifies the path or
#' paths that you are interested in
#' @param kegg_enzyme a kegg_tbl generated from the get_kegg_enzymes function. If
#' NULL (default) the kegg_enzymes dataset in this package will be used
#' @param kegg_orthology a kegg_tbl generated from the get_kegg_orthology function.
#' If NULL (default) the kegg_orthologies dataset in this package will be used.
#' @param kegg_module a kegg_tbl generated from the get_kegg_module function.
#' If NULL (default) the kegg_modules dataset in this package will be used
#' @param kegg_pathway a kegg_tbl generated from the get_kegg_pathway function.
#' If NULL (default) the kegg_pathways function will be used.
#' @param strict whether strict matching should be used when matching the pathway_name
#' provided to the pathways in the kegg_pathway dataset. If FALSE (default), regex will
#' be used to match pathways to the provided options. If TRUE, the provided pathways
#' will have to be the same as they are in the dataset (case insensitive)
#'
#' @export
keggtap <- function(pathway_name, kegg_enzyme = NULL,
                    kegg_orthology = NULL, kegg_module = NULL,
                    kegg_pathway = NULL, strict = FALSE){

  if (!match_strict){
    pathway_match <- paste(tolower(pathway), collapse = "|")
    pathways <- kegg_pathway[stringr::str_detect(tolower(kegg_pathway$pathway), pathway_match), ]
  } else {
    pathway_match <- tolower(pathway)
    pathways <- kegg_pathway[tolower(kegg_pathway$pathway) %in% pathway_match, ]
  }

  if (nrow(pathways) == 0){
    stop("There are no pathways that match your search", call. = FALSE)
  }

  pathway_enzymes <- get_pathway_enzymes(pathway_name, kegg_enzyme = kegg_enzyme, kegg_pathway = kegg_pathway, strict = strict)
  pathway_orthology <- get_pathway_orthologies(pathway_name, kegg_orthology = kegg_orthology, kegg_pathway = kegg_pathway, strict = strict)
  pathway_modules <- get_pathway_modules(pathway_name, kegg_module = kegg_module, kegg_pathway = kegg_pathway, strict = strict)

  out <- structure(
    list(
      pathway = pathway,
      enzyme = pathway_enzymes,
      orthology = pathway_orthology,
      module = pathway_modules
    ),
    class = "keggtap"
  )

  return(out)

}


#' @export
is_keggtap <- function(x){
  inherits(x, "keggtap")
}

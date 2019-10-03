# function to create keggtap object

#' Create a keggtap object
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




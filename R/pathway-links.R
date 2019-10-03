#' Query all modules associated with a given pathway
#'
#' @param pathway_name a character string specifying a certain pathway. Regular
#' expressions are used to match the pathway so full pathway names do not have
#' to be given
#'
#' @details the value passed to pathway_name is searched against the pathway
#' column in the kegg_pathways dataset
#'
#' @export
get_pathway_modules <- function(pathway_name, kegg_module = NULL, kegg_pathway = NULL, strict = FALSE) {
  if (is.null(kegg_pathway)){
    kegg_pathway <- KEGGerator::kegg_pathways
  } else {
    if (!is_kegg_tbl(kegg_pathway, "pathway")){
      stop("kegg_pathway must be a kegg_tbl with columns pathway and pathway_id", call. = FALSE)
    }
  }

  if (is.null(kegg_module)){
    kegg_module <- KEGGerator::kegg_modules
  } else {
    if (!is_kegg_tbl(kegg_module, "module")){
      stop("kegg_module must be a kegg_tbl with columns module and module_id", call. = FALSE)
    }
  }

  pathway <- kegg_pathway %>%
    filter_pathway(pathway_name, strict = strict) %>%
    # dplyr::filter(stringr::str_detect(pathway, pathway_name)) %>%
    dplyr::mutate(
      module = purrr::map(pathway_id, ~{
        mods <- KEGGREST::keggLink("module", .x)

        output <- tibble::tibble(
          module_id = mods
        ) %>%
          dplyr::left_join(kegg_module, by = "module_id")

        return(output)
      })
    ) %>%
    tidyr::unnest(module)

  return(pathway)
}

## THIS IS USED IN THE get_genome_enzymes AND get_genome_orthologies FUNCTIONS
## AND WILL NEED TO BE MODIFIED ONCE CLASSES ARE ADDED TO THESE OUTPUTS IF THAT
## ROUTE IS PURSUED


#' Query all enzymes associated with a given pathway
#'
#' @param pathway_name a character string specifying a certain pathway. Regular
#' expressions are used to match the pathway so full pathway names do not have
#' to be given.
#'
#' @details the value passed to pathway_name is searched against the pathway
#' column in the kegg_pathways dataset
#'
#' @export
get_pathway_enzymes <- function(pathway_name, kegg_enzyme = NULL, kegg_pathway = NULL, strict = FALSE) {
  if (is.null(kegg_pathway)){
    kegg_pathway <- KEGGerator::kegg_pathways
  } else {
    if (!is_kegg_tbl(kegg_pathway, "pathway")){
      stop("kegg_pathway must be a kegg_tbl with columns pathway and pathway_id", call. = FALSE)
    }
  }

  if (is.null(kegg_enzyme)){
    kegg_enzyme <- KEGGerator::kegg_enzymes
  } else if (!is_kegg_tbl(kegg_enzyme, "enzyme")){
    stop("kegg_enzyme must be a kegg_tbl with columns enzyme and enzyme_id", call. = FALSE)
  }

  pathway <- kegg_pathway %>%
    filter_pathway(pathway_name, strict = strict) %>%
    # dplyr::filter(stringr::str_detect(pathway, pathway_name)) %>%
    dplyr::mutate(
      enzyme = purrr::map(pathway_id, ~{
        enzymes <- KEGGREST::keggLink("enzyme", .x)

        output <- tibble::tibble(
          enzyme_id = enzymes
        ) %>%
          dplyr::left_join(kegg_enzyme, by = "enzyme_id")

        return(output)
      })
    ) %>%
    tidyr::unnest(enzyme)

  return(pathway)
}

#' Query all orthologies associated with a given pathway
#'
#' @param pathway_name a character string specifying a certain pathway. Regular
#' expressions are used to match the pathway so full pathway names do not have
#' to be given.
#'
#' @details the value passed to pathway_name is searched against the pathway
#' column in the kegg_pathways dataset
#'
#' @export
get_pathway_orthologies <- function(pathway_name, kegg_orthology = NULL, kegg_pathway = NULL, strict = FALSE){
  if (is.null(kegg_pathway)){
    kegg_pathway <- KEGGerator::kegg_pathways
  } else {
    if (!is_kegg_tbl(kegg_pathway, "pathway")){
      stop("kegg_pathway must be a kegg_tbl with columns pathway and pathway_id", call. = FALSE)
    }
  }

  if (is.null(kegg_orthology)){
    kegg_orthology <- KEGGerator::kegg_orthologies
  } else{
    if (!is_kegg_tbl(kegg_orthology, "orthology")){
      stop("kegg_orthology must be a kegg_tbl with columns orthology and orthology_id", call. = FALSE)
    }
  }


  pathway <- kegg_pathway %>%
    filter_pathway(pathway_name, strict = strict) %>%
    # dplyr::filter(stringr::str_detect(pathway, pathway_name)) %>%
    dplyr::mutate(
      orthology = purrr::map(pathway_id, ~{
        orth <- KEGGREST::keggLink("orthology", .x)

        output <- tibble::tibble(
          orthology_id = orth
        ) %>%
          dplyr::left_join(kegg_orthology, by = "orthology_id")

        return(output)
      })
    ) %>%
    tidyr::unnest(orthology)

  return(pathway)
}


filter_pathway <- function(kegg_pathway, pathway_name, strict){
  if (!is_kegg_tbl(kegg_pathway, "pathway")){
    stop("kegg_pathway must be a kegg_tbl with columns pathway and pathway_id", call. = FALSE)
  }

  if (strict){
    out <- dplyr::filter(kegg_pathway, pathway %in% pathway_name)
  } else {
    if (length(pathway_name) > 1) pathway_name <- paste(pathway_name, collapse = "|")
    out <- dplyr::filter(kegg_pathway, stringr::str_detect(pathway, pathway_name))
  }

  return(out)
}

#' Filter OTUs and kegg organisms
#'
#' @export
filter_orgs <- function(data, uncertainty, pathway_name, pathways, verbose, progress){
  UseMethod("filter_orgs")
}

#' @describeIn filter_orgs method for keggerator objects
#' @export
filter_orgs.keggerator <- function(data, uncertainty = 1, pathway_name = NULL,
                                   pathways = NULL, verbose = FALSE, progress = TRUE){

  if (is.null(data$orgs_tbl)) stop("orgs_tbl slot of keggerator object is NULL, have you run orgs_tibble()?", call. = FALSE)
  if(!is_orgs_tbl(data$orgs_tbl)) stop("orgs_tbl slow is not of class orgs_tbl. Did you generate it with orgs_tibble()?", call. = FALSE)
  if (is.null(data$orgs_id)) stop("orgs_id slot of keggerator object is NULL, have you run get_org_ids()?", call. = FALSE)
  if (!is_orgs_id(data$orgs_id)) stop("orgs_id slot is not of class orgs_id. Did you generate it with get_org_ids()?", call. = FALSE)

  filtered_orgs <- filter_orgs_internal(
    orgs_id = data$orgs_id,
    orgs_tbl = data$orgs_tbl,
    uncert_tbl = data$total_uncert,
    uncertainty = uncertainty,
    pathway_name = pathway_name,
    pathways = pathways,
    verbose = verbose,
    progress = progress
  )

  data$orgs_filt <- filtered_orgs

  return(data)

}

filter_orgs.orgs_list <- function(data, uncertainty = 1, pathway_name = NULL,
                                  pathways = NULL, verbose = FALSE, progress = TRUE){
  if (is.null(data$orgs_tbl)) stop("orgs_tbl slot of orgs_list object is NULL, have you run orgs_tibble()?", call. = FALSE)
  if(!is_orgs_tbl(data$orgs_tbl)) stop("orgs_tbl slow is not of class orgs_tbl. Did you generate it with orgs_tibble()?", call. = FALSE)
  if (is.null(data$orgs_id)) stop("orgs_id slot of orgs_list object is NULL, have you run get_org_ids()?", call. = FALSE)
  if (!is_orgs_id(data$orgs_id)) stop("orgs_id slot is not of class orgs_id. Did you generate it with get_org_ids()?", call. = FALSE)

  filtered_orgs <- filter_orgs_internal(
    orgs_id = data$orgs_id,
    orgs_tbl = data$orgs_tbl,
    uncert_tbl = data$total_uncert,
    uncertainty = uncertainty,
    pathway_name = pathway_name,
    pathways = pathways,
    verbose = verbose,
    progress = progress
  )

  data$orgs_filt <- filtered_orgs

  return(data)
}


filter_orgs_internal <- function(orgs_id, orgs_tbl, uncert_tbl = NULL, uncertainty = 1,
                                 pathway_name = NULL, pathways = NULL, verbose = FALSE,
                                 progress = TRUE){

  if (!is.null(uncert_tbl)){
    out <- filter_orgs_uncert(orgs_id = orgs_id, orgs_tbl = orgs_tbl,
                              uncert_tbl = uncert_tbl, uncertainty = uncertainty)

    n_removed_uncert <- nrow(orgs_id) - nrow(out)
    # ids_removed_uncert <- orgs_id$genome_id[!orgs_id$genome_id %in% out$genome_id]

    if (verbose | progress){
      cat(crayon::red(n_removed_uncert), " organisms removed because uncertainty was greater than ", crayon::red(uncertainty), "\n")
    }
  }





  if (!is.null(pathway_name)){
    n_current <- nrow(out)

    out <- filter_orgs_pathway(orgs_id = out, pathway_name = pathway_name, pathways = pathways)

    n_removed_path <- n_current - nrow(out)

    if (verbose | progress){
      cat(crayon::red(n_removed_path), " organisms removed because they were not linked with the ",
          crayon::red(pathway_name), " pathway in KEGG\n")
    }

  }

  if (is.null(uncert_tbl) & is.null(pathway_name)){
    warning("uncert_tbl and pathway_name were both NULL. No filtration was performed and the original object will be returned", call. = FALSE)
    return(orgs_id)
  } else {
    attr(out, "filtered") <- TRUE
    return(out)
  }

}

filter_orgs_uncert <- function(orgs_id, orgs_tbl, uncert_tbl, uncertainty = 1){

  if (!is_orgs_id(orgs_id)) stop("orgs_id must be of class orgs_id", call. = FALSE)
  if (!is_uncert_tbl(uncert_tbl)) stop("uncert_tbl must be of class uncert_tbl", call. = FALSE)

  otu_id_keep <- uncert_tbl %>%
    dplyr::filter(total_uncert <= uncertainty) %>%
    dplyr::pull(otu_id)

  otu_keep <- orgs_tbl %>%
    dplyr::filter(otu_id %in% otu_id_keep)

  org_keep <- otu_keep %>%
    dplyr::pull(genome)

  id_keep <- orgs_id %>%
    dplyr::filter(genome %in% org_keep)

  output <- id_keep

  attr(output, "genome_removed_uncert") <- orgs_id$genome_id[!orgs_id$genome %in% org_keep]

  return(output)


}



filter_orgs_pathway <- function(orgs_id, pathway_name, pathways = NULL, verbose = FALSE, progress = TRUE){
  if (is.null(pathways)){
    pathways <- kegg_pathways
  } else{
    if (!is_kegg_tbl(pathways, "pathway")){
      stop("pathways must be a kegg_tbl with the columns `pathway` and `pathway_id`", call. = FALSE)
    }
  }

  if (!any(stringr::str_detect(stringr::str_to_lower(pathways$pathway), stringr::str_to_lower(pathway_name)))){
    stop("pathway_name had no matches in pathways dataset", call. = FALSE)
  }

  pathway_ids <- pathways %>%
    dplyr::filter(stringr::str_detect(stringr::str_to_lower(pathway), stringr::str_to_lower(pathway_name))) %>%
    dplyr::mutate(pathway_id = stringr::str_replace(pathway_id, "path:map", "")) %>%
    dplyr::pull(pathway_id)

  if (progress){
    p <- dplyr::progress_estimated(nrow(orgs_id), 10)
  }


  output <- orgs_id %>%
    dplyr::mutate(
      genome_id = stringr::str_replace(genome_id, "gn:", ""),
      pathway = purrr::map2(genome_id, genome_desc, ~{

        paths = kegg_link_safe("pathway", .x)

        if (verbose){
          cat(
            "Linking pathways to ", crayon::red(.y), ": ",
            crayon::blue(length(paths)), " linked\n", sep = ""
          )
        }

        if (progress){
          p$tick()$print()
        }

        return(paths)
      }),
      # pathway = purrr::map(pathway, ~purrr::discard(.x, function(x){!stringr::str_detect(x, pathway_ids)}))
      pathway = purrr::map(pathway, ~{
        prefix <- stringr::str_replace_all(.x, "[0-9]", "") %>%
          unique()

        paths <- stringr::str_replace_all(.x, "[A-z]|\\:", "") %>%
          .[. %in% pathway_ids] %>%
          unique()

        output <- paste0(prefix, paths)
        return(output)
      })
    ) %>%
    dplyr::filter(purrr::map_lgl(pathway, ~{length(.x) > 0}))

  if (verbose & progress){
    p$stop()
  }

  attr(output, "genome_removed_pathway") <- orgs_id$genome_id[!orgs_id$stringr::str_replace(genome_id, "gn:", "") %in% output$genome_id]

  return(output)

}






#' Filter out genomes by pathway
#'
#' @param data a tibble that has genome_ids - likely the output from get_genome_id()
#' @param pathway_name the name of the pathway that you wish to filter
#' @param pathways a kegg_tbl dataset of available kegg pathways. This dataset
#' is produced from the get_kegg_pathway function. If it is NULL (default) the
#' KEGGerator::kegg_pathways dataset will be used
#'
#' @export
filter_genome_by_pathway <- function(data, pathway_name, pathways = NULL){
  if (is.null(pathways)){
    pathways <- kegg_pathways
  } else{
    if (!is_kegg_tbl(pathways, "pathway")){
      stop("pathways must be a kegg_tbl with the columns `pathway` and `pathway_id`", call. = FALSE)
    }
  }

  pathway_ids <- pathways %>%
    dplyr::filter(stringr::str_detect(pathway, pathway_name)) %>%
    dplyr::mutate(pathway_id = stringr::str_replace(pathway_id, "path:map", "")) %>%
    dplyr::pull(pathway_id)

  output <- data %>%
    dplyr::mutate(
      genome_id = stringr::str_replace(genome_id, "gn:", ""),
      pathway = purrr::map(genome_id, ~{
        kegg_link_safe("pathway", .x)
      }),
      # pathway = purrr::map(pathway, ~purrr::discard(.x, function(x){!stringr::str_detect(x, pathway_ids)}))
      pathway = purrr::map(pathway, ~{
        prefix <- stringr::str_replace_all(.x, "[0-9]", "") %>%
          unique()

        paths <- stringr::str_replace_all(.x, "[A-z]|\\:", "") %>%
          .[. %in% pathway_ids] %>%
          unique()

        output <- paste0(prefix, paths)
        return(output)
      })
    ) %>%
    dplyr::filter(purrr::map_lgl(pathway, ~{length(.x) > 0}))

  return(output)
}

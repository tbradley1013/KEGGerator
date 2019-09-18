#' Calculate total uncertainty
#'
#' @export
keggerator_uncertainty <- function(data){
  UseMethod("keggerator_uncertainty")
}

#' @describeIn keggerator_uncertainty method for keggerator objects
#' @export
keggerator_uncertainty.keggerator <- function(data){
  spec <- data$species_uncert
  kegg <- data$kegg_uncert

  if (is.null(spec)) stop("the `species_uncert` slot is empty in your keggerator object, have you run orgs_tibble() yet?", call. = FALSE)
  if (is.null(kegg)) stop("the `kegg_uncert` slot is empty in your keggerator object, have you run get_org_ids() yet?", call. = FALSE)

  total_uncert <- uncertainty_internal(spec, kegg)

  data$total_uncert <- total_uncert

  return(data)
}

#' @describeIn keggerator_uncertainty method for orgs_list objects
#' @export
keggerator_uncertainty.orgs_list <- function(data){
  spec <- data$species_uncert
  kegg <- data$kegg_uncert

  if (is.null(spec)) stop("the `species_uncert` slot is empty in your orgs_list object, have you run orgs_tibble() yet?", call. = FALSE)
  if (is.null(kegg)) stop("the `kegg_uncert` slot is empty in your orgs_list object, have you run get_org_ids() yet?", call. = FALSE)

  data$total_uncert <- total_uncert
  return(data)
}



uncertainty_internal <- function(spec_uncert, kegg_uncert){

    kegg_uncert_otu <- kegg_uncert %>%
    dplyr::group_by(otu_id) %>%
    dplyr::summarize(total_kegg_hits = sum(n_hits)) %>%
    dplyr::ungroup()

  output <- spec_uncert %>%
    dplyr::left_join(kegg_uncert_otu, by = "otu_id") %>%
    dplyr::mutate(total_uncert = dplyr::if_else(total_kegg_hits == 0, 1, 1-(1/max(total_kegg_hits, n_spec)))) %>%
    dplyr::select(otu_id, n_spec, total_kegg_hits, total_uncert)

  return(output)
}

#' Calculate total uncertainty
#'
#' @export
keggerator_uncertainty <- function(data){
  UseMethod("keggerator_uncertainty")
}


keggerator_uncertainty.keggerator <- function(data){
  spec <- data$species_uncert
  kegg <- data$kegg_uncert

  if (is.null(spec)) stop("the `species_uncert` slot is empty in your keggerator object, have you run orgs_tibble() yet?", call. = FALSE)
  if (is.null(kegg)) stop("the `kegg_uncert` slot is empty in your keggerator object, have you run get_org_ids() yet?", call. = FALSE)
}



uncertainty_internal <- function(spec_uncert, kegg_uncert){

    kegg_uncert_otu <- kegg_uncert %>%
    dplyr::group_by(otu_id) %>%
    dplyr::summarize(total_kegg_hits = sum(n_hits)) %>%
    dplyr::ungroup()

  output <- spec_uncert %>%
    dplyr::left_join(kegg_uncert_otu, by = "otu_id") %>%
    dplyr::mutate(uncert = dplyr::if_else(total_kegg_hits == 0, 1, 1-(1/max(total_kegg_hits, n_spec))))

  return(output)
}

#' Query all modules related to each enzyme in a dataset
#'
#' @param data a tibble that is output from get_pathway_enzymes()
#'
#' @export
get_enzyme_modules <- function(data, kegg_module = NULL){
  if (!"tbl_df" %in% class(data)) stop("data must be of class tbl_df")

  if (!"enzyme_id" %in% colnames(data)) stop("data must contian a column named 'enzyme_id'")

  if (is.null(kegg_module)){
    kegg_module <- KEGGerator::kegg_modules
  } else {
    if (!is_kegg_tbl(kegg_module, "module")){
      stop("kegg_module must be a kegg_tbl with columns module and module_id", call. = FALSE)
    }
  }

  output <- data %>%
    dplyr::mutate(
      module = purrr::map(enzyme_id, ~{
        mods <- KEGGREST::keggLink("module", .x)

        if (length(mods) > 0) {
          output <- tibble::tibble(
            module_id = mods
          ) %>%
            dplyr::left_join(kegg_module, by = "module_id")
        } else {
          output <- tibble::tibble(
            module_id = "None",
            module = "None"
          )
        }

        return(output)
      })
    ) %>%
    tidyr::unnest(module)

  return(output)
}


#' Query all modules related to each orthology in a dataset
#'
#' @param data a tabble that is output from get_pathway_orthologies()
#'
#' @export
get_orthology_modules <- function(data){
  if (!"tbl_df" %in% class(data)) stop("data must be of class tbl_df")

  if (!"orthology_id" %in% colnames(data)) stop("data must contain a column named 'orthology_id'")

  output <- data %>%
    dplyr::mutate(
      module = purrr::map(orthology_id, ~{
        mods <- KEGGREST::keggLink("module", .x)

        if (length(mods) > 0){
          output <- tibble::tibble(
            module_id = mods
          ) %>%
            dplyr::left_join(kegg_modules, by = "module_id")
        } else {
          output <- tibble::tibble(
            module_id = "None",
            module = "None"
          )
        }

        return(output)
      })
    ) %>%
    tidyr::unnest(module)

  return(output)
}

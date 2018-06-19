#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

# safe keggGet function
kegg_get_safe <- purrr::possibly(KEGGREST::keggGet, otherwise = list(NA_character_))

# safe keggLink function
kegg_link_safe <- purrr::possibly(KEGGREST::keggLink, otherwise = NA_character_)

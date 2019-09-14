#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

# safe keggGet function
kegg_get_safe <- purrr::possibly(KEGGREST::keggGet, otherwise = list(NA_character_))

# safe keggLink function
kegg_link_safe <- purrr::possibly(KEGGREST::keggLink, otherwise = NA_character_)

# safe keggList function
kegg_list_safe <- purrr::possibly(KEGGREST::keggList, otherwise = NA_character_)

# safe keggFind function
kegg_find_safe <- purrr::possibly(KEGGREST::keggFind, otherwise = NA_character_)

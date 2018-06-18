kegg_get_safe <- purrr::possibly(KEGGREST::keggGet, otherwise = list(NA_character_))

kegg_link_safe <- purrr::possibly(KEGGREST::keggLink, otherwise = NA_character_)

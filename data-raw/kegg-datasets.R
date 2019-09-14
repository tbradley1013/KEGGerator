library(KEGGREST)
library(tibble)


kegg_pathways_list <- KEGGREST::keggList("pathway")

kegg_pathways <- tibble::tibble(
  pathway_id = names(kegg_pathways_list),
  pathway = kegg_pathways_list
)

class(kegg_pathways) <- c("kegg_tbl", class(kegg_pathways))

kegg_module_list <- KEGGREST::keggList("module")

kegg_modules <- tibble::tibble(
  module_id = names(kegg_module_list),
  module = kegg_module_list
)

class(kegg_modules) <- c("kegg_tbl", class(kegg_modules))

kegg_enzyme_list <- KEGGREST::keggList("enzyme")

kegg_enzymes <- tibble::tibble(
  enzyme_id = names(kegg_enzyme_list),
  enzyme = kegg_enzyme_list
)

class(kegg_enzymes) <- c("kegg_tbl", class(kegg_enzymes))

kegg_orthology_list <- KEGGREST::keggList("orthology")

kegg_orthologies <- tibble::tibble(
  orthology_id = names(kegg_orthology_list),
  orthology = kegg_orthology_list
)

class(kegg_orthologies) <- c("kegg_tbl", class(kegg_orthologies))


usethis::use_data(kegg_pathways, overwrite = TRUE, compress = "gzip")
usethis::use_data(kegg_modules, overwrite = TRUE, compress = "gzip")
usethis::use_data(kegg_enzymes, overwrite = TRUE, compress = "gzip")
usethis::use_data(kegg_orthologies, overwrite = TRUE, compress = "gzip")

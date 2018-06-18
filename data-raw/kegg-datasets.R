library(KEGGREST)
library(tibble)


kegg_pathways_list <- KEGGREST::keggList("pathway")

kegg_pathways <- tibble::tibble(
  pathway_id = names(kegg_pathways_list),
  pathway = kegg_pathways_list
)


kegg_module_list <- KEGGREST::keggList("module")

kegg_modules <- tibble::tibble(
  module_id = names(kegg_module_list),
  module = kegg_module_list
)

kegg_enzyme_list <- KEGGREST::keggList("enzyme")

kegg_enzymes <- tibble::tibble(
  enzyme_id = names(kegg_enzyme_list),
  enzyme = kegg_enzyme_list
)

kegg_orthology_list <- KEGGREST::keggList("orthology")

kegg_orthologies <- tibble::tibble(
  orthology_id = names(kegg_orthology_list),
  orthology = kegg_orthology_list
)


devtools::use_data(kegg_pathways, overwrite = TRUE, compress = "gzip")
devtools::use_data(kegg_modules, overwrite = TRUE, compress = "gzip")
devtools::use_data(kegg_enzymes, overwrite = TRUE, compress = "gzip")
devtools::use_data(kegg_orthologies, overwrite = TRUE, compress = "gzip")

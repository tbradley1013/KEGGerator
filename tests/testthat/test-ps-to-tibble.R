library(phyloseq)

ps_file <- system.file("data", "GlobalPatterns.Rdata", package = "phyloseq")
ps <- filter_taxa(load(ps_file), function(x)sum(x)>10000, prune = TRUE)

test_that("conversion from phyloseq to keggerator components works", {

  ps_otu <- otu_tibble(ps)
  ps_otu_ref <- otu_ref(ps)
  ps_tax <- tax_tibble(ps)



  # testing that the otu ids work correctly
  expect_identical(ps_otu_ref, dplyr::select(ps_otu, otu_id, otu))
  expect_identical(ps_otu_ref, dplyr::select(ps_tax, otu_id, otu))
})

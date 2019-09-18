#===============================================================================
# workflow to help test interactively during development.... need to converrt
# to unit tests
#
# Tyler Bradley
# 2019-09-16
#===============================================================================

library(phyloseq)
library(KEGGerator)
library(tidyverse)

ps <- read_rds("../WWTP_Impact_on_Stream/data/ps_orig.RDS")

# ps <- subset_taxa(ps, !is.na(Species))

ps <- filter_taxa(ps, function(x)sum(x)>1000, prune = TRUE)

# gerenate a
ps_kegg <- as_keggerator(ps)

ps_kegg <- orgs_tibble(ps_kegg)

# orgs_id <- get_org_ids(ps_kegg$orgs_tbl)

# get org ids inside keggerator object
ps_kegg <- get_org_ids(ps_kegg)


ps_genomes <- genomes_tibble(ps)


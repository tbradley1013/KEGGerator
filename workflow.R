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

# ps <- readr::read_rds("../WWTP_Impact_on_Stream/data/ps_orig.RDS")
#
# ps <- phyloseq::subset_taxa(ps, !is.na(Species))
#
# ps <- phyloseq::filter_taxa(ps, function(x)sum(x)>500, prune = TRUE)
#
# # gerenate a
# ps_kegg <- as_keggerator(ps)
#
# ps_kegg <- orgs_tibble(ps_kegg)
#
# # orgs_id <- get_org_ids(ps_kegg$orgs_tbl)
#
# # get org ids inside keggerator object
# ps_kegg <- get_org_ids(ps_kegg)
#
# # ps_kegg <- filter_orgs()
#
# ps_kegg <- get_enzymes(ps_kegg)
#
# ps_kegg <- get_orthologies(ps_kegg)
#
# write_rds(ps_kegg, "workflow-keggerator.rds", compress = "xz", compression = 9L)

ps_kegg <- read_rds("workflow-keggerator.rds")

nit_enzymes <- get_pathway_enzymes("Nitrogen metabolism")

nit_orths <- get_pathway_orthologies("Nitrogen metabolism")
nit_mods <- get_pathway_modules("Nitrogen metabolism")



pan_ko <- ps_kegg$orgs_tbl %>%
  left_join(ps_kegg$orgs_orthologies, by = "genome") %>%
  filter(!is.na(genome_id)) %>%
  group_by(otu_id, orthology_id) %>%
  summarize(n_ko = n()) %>%
  ungroup() %>%
  left_join(ps_kegg$total_uncert %>%
              select(otu_id, total_kegg_hits),
            by = "otu_id") %>%
  mutate(
    pan = n_ko/total_kegg_hits
  )

orgs_prob <- ps_kegg$total_uncert %>%
  mutate(prob = ifelse(n_spec >= total_kegg_hits, 1/n_spec, 1/total_kegg_hits)) %>%
  select(otu_id, prob) %>%
  left_join(pan_ko, by = "otu_id") %>%
  filter(!is.na(pan), !is.na(prob)) %>%
  mutate(overall_prob = prob*pan) %>%
  select(otu_id, orthology_id, overall_prob)

ps_kegg$otu_tbl %>% mutate_at(vars(-c(otu_id, otu)), list(~./sum(.)))



# ps_kegg$orgs_tbl %>%
#   left_join(ps_kegg$otu_tbl, by = "otu_id") %>%
#   gather(key = sample, value = count, -c(otu_id, genome, otu)) %>%
#   group_by(genome) %>%
#   mutate(count = count/median(count)) %>%
#   ungroup() %>%
#   group_by(sample) %>%
#   mutate(count = count/sum(count)) %>%
#   ungroup() %>%
#   spread(key = sample, value = count) %>%
#   select(-genome) %>%
#   distinct()



ps_kegg$orgs_tbl %>%
  left_join(ps_kegg$otu_tbl, by = "otu_id") %>%
  mutate_at(vars(-c(otu_id, otu, genome)), list(~./sum(.))) %>%
  inner_join(orgs_prob, by = "otu_id") %>%
  mutate_at(vars(-c(otu_id, genome, otu, orthology_id, overall_prob)), list(~.*overall_prob)) %>%
  group_by(orthology_id) %>%
  summarize_at(vars(-c(otu_id, genome, otu, overall_prob)), sum)

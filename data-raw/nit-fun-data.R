#===============================================================================
# This will include dataset of all of the functionalities of the Nitrogen
# metabolism and what functionalities are required for each individual step
# in the nitrogen cycle.
#===============================================================================

# This will read in a dataset that was hand editted in excel to group each of
# the individual orthologies. Some orthologies are in the dataset twice if they
# serve more than one purpose.
nit_orth_funcs <- read_csv("data-raw/nit-orth-funcs.csv")

# This tibble puts the function ids (created by me) in a tibble with their
# long hand meaning (again written by me)
function_key <- tibble::tribble(
  ~func_id,             ~func_meaning,
  "NO3R",       "Nitrate reduction",
  "NO2R",       "Nitrite reduction",
  "NO2O",       "Nitrite oxidation",
  "N2OR", "Nitrous oxide reduction",
  "NF",       "Nitrogen fixation",
  "NOR",  "Nitric oxide reduction",
  "NO3A",    "Nitrate assimilation",
  "NH2OHO", "Hydroxylamine oxidation",
  "NH3O",       "Ammonia oxidation",
  "NH3A1",    "Ammonia assimilation to glutamate",
  "GLTS",                 "Glutamate synthesis",
  "NH3A2",   "Ammonia assimilation to glutamine",
  "HYDR",             "Hydrazine dehydrogenase",
  "NOR2", "Nitric oxide reduction to hyrdazine",
  "NH3O2",      "Ammonia oxidation to hydrazine"
)

nitrogen_metabolism_orth_fun <- nit_orth_funcs %>%
  dplyr::rename(orthology_id = orthology) %>%
  dplyr::left_join(function_key, by = "func_id")

# This tibble shows each of the modules (not just the ones from KEGG) for the
# Nitrogen metabolism pathway and all of the functions (as defined aboved) that
# are required for the each one
nitrogen_metabolism_fun_groups <- tribble(
  ~modules_name, ~funcs_in_module,
  "Assimilatory nitrate reduction, nitrate => ammonia", list("NO3R", "NO2R"),
  "Dissimilatory nitrate reduction, nitrate => ammonia", list("NO3R", "NO2R"),
  "Denitrification, nitrate => nitrogen", list("NO3R", "NO2R", "NOR", "N2OR"),
  "Complete nitrification, comammox, ammonia => nitrite => nitrate", list("NH3O", "NH2OHO", "NO2O"),
  "Nitrogen fixation, nitrogen => ammonia", list("NF"),
  "Nitrate assimilation", list("NO3A"),
  "Nitrification, ammonia => nitrite", list("NH3O", "NH2OHO"),
  "Nitrate/nitrite transport system", list("NO3A"),
  "Ammonia assimilation, ammonia => glutamate", list("NH3A1"),
  "Ammonia assimilation, ammonia => glutamine => glutamate", list("NH3A2", "GLTS"),
  "Ammamox", list("HYDR", "NOR2", "NH3O2")

)


devtools::use_data(nitrogen_metabolism_orth_fun, overwrite = TRUE)
devtools::use_data(nitrogen_metabolism_fun_groups)

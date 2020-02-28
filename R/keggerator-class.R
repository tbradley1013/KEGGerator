#' A keggerator object
#'
#' A class of object that contains all of the
#' required components for analysis with
#' KEGGerator
#'
#' @return an object of class keggerator
#'
#' @export
keggerator <- R6::R6Class(
  "keggerator",
  public = list(
    ps = NULL,
    tax_tbl = NULL,
    otu_tbl = NULL,
    sam_tbl = NULL,
    otu_ref = NULL,
    orgs_tbl = NULL,
    species_uncert = NULL,
    orgs_id = NULL,
    kegg_uncert = NULL,
    total_uncert = NULL,
    orgs_filt = NULL,
    orgs_enzymes = NULL,
    orgs_orthologies = NULL,

    initialize = function(ps, .defaults = TRUE){
      stopifnot(inherits(ps, "phyloseq"))

      cat("Intializing KEGGerator...\n")
      tictoc::tic("Finished Converting phyloseq object")
      self$ps <- ps
      self$tax_tbl <- verbose("Converting Tax table", tax_tibble, self$ps)
      self$otu_tbl <- verbose("Converting OTU table", otu_tibble, self$ps)
      self$sam_tbl <- verbose("Converting Sample table", sam_tibble, self$ps)
      self$otu_ref <- verbose("Creating OTU reference table", otu_ref, self$ps)
      tictoc::toc()

      if (.defaults){
        verbose("Generating organisms table", self$make_orgs_tbl)
      }

      invisible(self)
    },

    #' @description
    #' Extract organisms from tax table or tax tibble
    #'
    #' @param data either a phyloseq object or the output of tibble_tax
    #' @param drop_taxa logical; should the taxonomy be removed from the output
    #' @param strict logical; if set to TRUE, it will not keep any taxa that have
    #' not been classified to the species level.
    #' @param sep what is used to separate multiple species in the species column.
    #' Must be valid regex targeting the seperator. Default is "\\/" to match the
    #' default separator used during taxonomic classification by \code{\link[dada2]{addSpecies}}
    #'
    #' @details
    #' When an object of class tax_tbl is passed to orgs_tibble, a list with two
    #' tibbles will be returned. The first tibble is the orgs tibble and the second
    #' is the uncertainty tibble which is comprised of the percent uncertainty [0-1]
    #' of the species level assignment of each otu. If the otu was not assigned to
    #' the species level than the uncertainty level is 1. If the otu was assigned to
    #' only a single species than the uncertainty level is 0. If the species level
    #' was assigned to N possible species, than the uncertainty level is 1-(1/N).
    #'
    #' If an object of class keggerator is passed to orgs_tibble than the same two
    #' tbls are returned, but rather than being in a list alone, they are added to
    #' the keggerator object that is given in the orgs_tbl and species_uncert
    #' slots, respectively.
    make_orgs_tbl = function(drop_taxa = TRUE, strict = FALSE, sep = "\\/"){
      data <- self$tax_tbl

      if (is.null(data)){
        stop("Tax table is not defined, was keggerator object initiated successfully?", call. = FALSE)
      }

      if (!ids_match(data)){
        warning("tax_tbl provided has not had otu_ids verified against otu_tbl and may cause issues later if they do not match.",
                "Run KEGGerator:::check_otu_id() on your phyloseq data to ensure OTUs match.", call. = FALSE)
      }

      # removing all otus that are not defined to at least the genus level
      taxa <- data %>%
        dplyr::mutate_if(is.factor, as.character) %>%
        dplyr::filter(!is.na(Genus)) %>%
        dplyr::mutate(Species = dplyr::if_else(is.na(Species), "", Species))

      # getting all otu_ids that are not defined at the species level
      spec_undef <- taxa$otu_id[taxa$Species == ""]

      # removing the taxa that do not have any species assigned
      if (strict){
        taxa <- dplyr::filter(taxa, Species != "")
      }

      # finding the max uncertainty in any of the species columns
      max_uncertainty <- taxa$Species %>%
        stringr::str_count(sep) %>%
        max()

      taxa_cols <- purrr::map_chr(1:(max_uncertainty + 1), ~paste0("genome", .x))

      orgs <- taxa %>%
        {suppressWarnings(tidyr::separate(., Species, into = taxa_cols, sep = sep))} %>%
        tidyr::gather(key = key, value = Species, tidyselect::all_of(taxa_cols)) %>%
        dplyr::filter(!is.na(Species)) %>%
        tidyr::unite("genome", Genus, Species, sep = " ") %>%
        dplyr::mutate(genome = stringr::str_trim(genome))

      if (drop_taxa) orgs <- dplyr::select(orgs, otu_id, genome)

      otu_species_uncert <- orgs %>%
        dplyr::group_by(otu_id) %>%
        dplyr::summarize(n_spec = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          n_spec = dplyr::if_else(otu_id %in% spec_undef, 0L, n_spec),
          spec_uncert = dplyr::if_else(otu_id %in% spec_undef, 1, 1-(1/n_spec))
        )

      class(orgs) <- c("orgs_tbl", class(orgs))
      class(otu_species_uncert) <- c("uncert_tbl", class(otu_species_uncert))

      self$orgs_tbl <- orgs
      self$species_uncert <- otu_species_uncert

      invisible(self)
    },

    #' @description
    #' Query KEGG genome id for all organims
    #'
    #' @param data dataset of genome names to get id for
    #'
    #' @details
    #' This took ~4.4 minutes to run for 400 organisms and will likely
    #' vary from system to system
    get_orgs_ids = function(verbose = TRUE, progress = TRUE){
      data <- self$orgs_tbl

      if (is.null(data)){
        warning("keggerator object does not have `orgs_tbl` defined yet, trying to define it now")
        self$make_orgs_tbl()
      }

      if (!is_orgs_tbl(data)) stop("data must have column named genome", call. = FALSE)

      unq_orgs <- dplyr::distinct(data, genome)

      if (progress){
        p <- dplyr::progress_estimated(nrow(unq_orgs), 10)
      }

      org_hits <- unq_orgs %>%
        dplyr::mutate(genome_query = purrr::map(genome, ~{
          query <- kegg_find_safe("genome", .x)

          if (verbose){
            cat(
              "Finding genomes that match ", crayon::red(.x), ": ",
              crayon::blue(length(query)), "\n", sep = ""
            )
          }

          if (progress){
            p$tick()$print()
          }

          return(query)
        })) %>%
        dplyr::filter(purrr::map(genome_query, length) > 0,
                      !is.na(genome_query)) %>%
        dplyr::mutate(
          genome_id = purrr::map(genome_query, ~names(.x))
        ) %>%
        tidyr::unnest(genome_query, genome_id) %>%
        dplyr::distinct() %>%
        tidyr::separate(genome_query, into = c("genome_name", "genome_desc"), sep = "; ")


      kegg_uncert <- org_hits %>%
        dplyr::group_by(genome) %>%
        dplyr::summarise(
          n_hits = sum(!is.na(genome_id))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          uncert = dplyr::if_else(n_hits == 0, 1, 1-(1/n_hits))
        ) %>%
        dplyr::inner_join(ps_kegg$orgs_tbl, by = "genome") %>%
        dplyr::arrange(otu_id)

      class(org_hits) <- c("orgs_id", class(org_hits))
      attr(org_hits, "filtered") <- FALSE
      class(kegg_uncert) <- c("uncert_tbl", class(kegg_uncert))

      self$orgs_id <- org_hits
      self$kegg_uncert <- kegg_uncert

      if (verbose & progress){
        p$stop()
      }

      invisible(self)

    },

    #' @description
    #' Query enzymes associated with organisms from KEGG
    #'
    #' @param pathway_enzymes a vector of pathway enzymes
    #' @param kegg_enzymes a kegg_tbl that has the columns enzyme and enzyme_id. This
    #' can be generated using the get_kegg_enzyme() function. If NULL (default) than
    #' the KEGGerator::kegg_enzymes dataset will be used
    #' @param verbose logical; if TRUE the number of enzyme links for each genome
    #' id will be shown as they are processed
    #' @param progress logical; if TRUE (default) than a progress bar will appear if
    #' the query takes longer than 10 seconds (it likely will if your data has more)
    #' than only a few genomes
    #'
    #' @details if `pathway_enzymes` is NULL (default) than all enzymes for each genome
    #' will be returned. If a vector of enzyme_ids is provided than only the enzymes related
    #' to each genome that are in this vector will be returned. This is useful if used in
    #' coordination with the output of `get_pathway_enzymes` to only include the enzymes
    #' that are related with the pathway of interest.
    get_enzymes = function(pathway_enzymes = NULL, kegg_enzymes = NULL, verbose = FALSE, progress = TRUE){
      orgs_id <- self$orgs_id

      if (is.null(orgs_id)){
        warning("orgs_id is not defined in the keggerator object, trying to define it now", call. = FALSE)
        self$get_orgs_id()
      }

      if (is.null(kegg_enzymes)){
        kegg_enzymes <- KEGGerator::kegg_enzymes
      } else{
        if (!is_kegg_tbl(kegg_enzymes, "enzyme")){
          stop("kegg_enzyme must be a kegg_tbl with columns enzyme and enzyme_id", call. = FALSE)
        }
      }

      if (!is_filtered(orgs_id)){
        warning("The orgs_id object provided has not been filtered, did you forget to run filter_orgs()?", call. = FALSE)
      }

      if (is.null(pathway_enzymes)){
        path_enzymes <- NULL
      } else if (is_keggtap(pathway_enzymes)){
        path_enzymes <- pathway_enzymes$enzyme$enzyme_id
      } else if (tibble::is_tibble(pathway_enzymes) | is.data.frame(pathway_enzymes)){
        if (!"enzyme_id" %in% colnames(pathway_enzymes)){
          stop("if pathway_enzymes is a tbl_df than it must have a column named enzyme_id")
        }

        path_enzymes <- pathway_enzymes$enzyme_id
      } else {
        path_enzymes <- pathway_enzymes
      }


      if (progress){
        p <- dplyr::progress_estimated(nrow(orgs_id), 10)
      }

      output <- orgs_id %>%
        dplyr::mutate(
          enzyme_id = purrr::map(genome_id, ~{
            id <- stringr::str_replace(.x, "gn:", "")

            enzymes <- kegg_link_safe("enzyme", id) %>%
              unique()

            n_hits <- length(enzymes)

            if (!is.null(path_enzymes)) {
              enzymes <- enzymes[enzymes %in% path_enzymes]
            }

            n_remain <- length(enzymes)

            if (verbose){
              if (!is.null(path_enzymes)){
                cat("Linked ", crayon::red(n_hits), " (", crayon::red(n_remain),  " within the pathway specified) enzymes to genome: ", crayon::blue(.x), "\n")
              } else {
                cat("Linked ", crayon::red(n_hits), " enzymes to genome: ", crayon::blue(.x), "\n")
              }

            }

            if (progress){
              p$tick()$print()
            }

            return(enzymes)
          })
        ) %>%
        tidyr::unnest(enzyme_id) %>%
        dplyr::left_join(kegg_enzymes, by = "enzyme_id")


      class(output) <- c("kegg_fun", class(output))
      attr(output, "query") <- "enzymes"


      if (progress & verbose){
        p$stop()
      }

      self$orgs_enzymes <- output
      invisible(self)
    },

    #' @description
    #' Query all orthologies related to each genome in the dataset
    #'
    #' @param data a tibble with genomes and genome_id or the output from `get_genome_id()`
    #' @param pathway_orthologies a vector of pathway orthologies
    #' @param kegg_orthology a kegg_tbl with the columns orthology and orthology_id
    #' which can be generated with the get_kegg_orthology function. If NULL (default)
    #' the KEGGerator::kegg_orthologies dataset will be used
    #' @param verbose logical; if TRUE the number of enzyme links for each genome
    #' id will be shown as they are processed
    #' @param progress logical; if TRUE (default) than a progress bar will appear if
    #' the query takes longer than 10 seconds (it likely will if your data has more)
    #' than only a few genomes
    #'
    #' @details if `pathway_orthologies` is NULL (default) than all orthologies for each genome
    #' will be returned. If a vector of orthology_ids is provided than only the orthologies related
    #' to each genome that are in this vector will be returned. This is useful if used in
    #' coordination with the output of `get_pathway_orthologies` to only include the orthologies
    #' that are related with the pathway of interest.
    get_orthologies = function(pathway_orthologies = NULL, kegg_orthology = NULL,
                               verbose = FALSE, progress = TRUE){
      orgs_id <- self$orgs_id

      if (is.null(orgs_id)){
        warning("orgs_id is not defined in the keggerator object, trying to define it now", call. = FALSE)
        self$get_orgs_id()
      }


      if (is.null(kegg_orthology)){
        kegg_orthology <- KEGGerator::kegg_orthologies
      } else{
        if (!is_kegg_tbl(kegg_orthology, "orthology")){
          stop("kegg_orthology must be a kegg_tbl with columns orthology and orthology_id", call. = FALSE)
        }
      }

      if (!is_filtered(orgs_id)){
        warning("The orgs_id object provided has not been filtered, did you forget to run filter_orgs()?", call. = FALSE)
      }

      if (is.null(pathway_orthologies)){
        path_orthology <- NULL
      } else if (is_keggtap(pathway_orthologies)){
        path_orthology <- pathway_orthologies$orthology$orthology_id
      } else if (tibble::is_tibble(pathway_orthologies) | is.data.frame(pathway_orthologies)){
        if (!"orthology_id" %in% colnames(pathway_orthologies)){
          stop("if pathway_orthologies is a tbl_df than it must have a column named enzyme_id")
        }

        path_orthology <- pathway_orthologies$enzyme_id
      } else {
        path_orthology <- pathway_orthologies
      }

      if (progress){
        p <- dplyr::progress_estimated(nrow(orgs_id), 10)
      }

      output <- orgs_id %>%
        dplyr::mutate(
          orthology_id = purrr::map(genome_id, ~{
            id <- stringr::str_replace(.x, "gn:", "")

            orths <- kegg_link_safe("orthology", id) %>%
              unique()

            n_hits <- length(orths)

            if (!is.null(path_orthology)) {
              orths <- orths[orths %in% path_orthology]
            }

            n_remain <- length(orths)

            if (verbose){
              if (!is.null(path_orthology)){
                cat("Linked ", crayon::red(n_hits), " (", crayon::red(n_remain),  " within the pathway specified) orthologies to genome: ", crayon::blue(.x), "\n")
              } else {
                cat("Linked ", crayon::red(n_hits), " orthologies to genome: ", crayon::blue(.x), "\n")
              }
            }

            if (progress){
              p$tick()$print()
            }

            return(orths)
          })
        ) %>%
        tidyr::unnest(orthology_id) %>%
        dplyr::left_join(kegg_orthology, by = "orthology_id")

      class(output) <- c("kegg_fun", class(output))
      attr(output, "query") <- "orthologies"

      if (progress & verbose){
        p$stop()
      }

      self$orgs_orthologies <- output
      invisible(self)
    }


  )
)

null_check_req <- function(x, .f){
  if (is.null(x)) return(FALSE)

  .f(x)
}

null_check_opt <- function(x, .f) {
  if (is.null(x)) return(TRUE)

  .f(x)
}


verbose <- function(msg, .f, ...){
  cat(glue::glue("{msg}"), "\n")
  tictoc::tic(glue::glue("Finished {msg}"))
  out <- .f(...)
  tictoc::toc()
  return(out)
}

#' Convert an object into class keggerator
#'
#' @param ps a phyloseq object to be converted to a keggerator object
#' @param ... additional arguements to be passed to \code{\link[KEGGerator]{keggerator}}
#'
#' @export
as_keggerator <- function(ps, ...){
  UseMethod("as_keggerator")
}

#' @describeIn as_keggerator method for phyloseq objects
#' @export
as_keggerator.phyloseq <- function(ps, ...){
  output <- keggerator$new(ps)
  return(output)
}




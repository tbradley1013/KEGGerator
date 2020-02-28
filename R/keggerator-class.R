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

# verbose <- function(msg, expr){
#   cat(glue::glue("{msg}"), "\n")
#   tictoc::tic(glue::glue("Finished {msg}"))
#   fn <- exprToFunction(expr)
#   out <- fn()
#   tictoc::toc()
#   return(out)
# }

verbose <- function(msg, .f, ...){
  cat(glue::glue("{msg}"), "\n")
  tictoc::tic(glue::glue("Finished {msg}"))
  out <- .f(...)
  tictoc::toc()
  return(out)
}

# # From Shiny
# makeFunction <- function(args = pairlist(), body, env = parent.frame()) {
#   eval(call("function", args, body), env)
# }
#
# exprToFunction <- function(expr, env=parent.frame(), quoted=FALSE) {
#   if (!quoted) {
#     expr <- eval(substitute(substitute(expr)), parent.frame())
#   }
#
#   # expr is a quoted expression
#   makeFunction(body=expr, env=env)
# }


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




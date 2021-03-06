% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orgs-tbl.R
\name{orgs_tibble}
\alias{orgs_tibble}
\alias{orgs_tibble.tax_tbl}
\alias{orgs_tibble.keggerator}
\alias{orgs_tibble.phyloseq}
\title{extract organisms from tax table or tax tibble}
\usage{
orgs_tibble(data, drop_taxa, strict, sep)

\method{orgs_tibble}{tax_tbl}(data, drop_taxa = TRUE, strict = FALSE,
  sep = "\\\\/")

\method{orgs_tibble}{keggerator}(data, drop_taxa = TRUE,
  strict = FALSE, sep = "\\\\/")

\method{orgs_tibble}{phyloseq}(data, drop_taxa = TRUE, strict = FALSE,
  sep = "\\\\/")
}
\arguments{
\item{data}{either a phyloseq object or the output of tibble_tax}

\item{drop_taxa}{logical; should the taxonomy be removed from the output}

\item{strict}{logical; if set to TRUE, it will not keep any taxa that have
not been classified to the species level.}

\item{sep}{what is used to separate multiple species in the species column.
Must be valid regex targeting the seperator. Default is "\\/" to match the
default separator used during taxonomic classification by \code{\link[dada2]{addSpecies}}}
}
\description{
extract organisms from tax table or tax tibble
}
\details{
When an object of class tax_tbl is passed to orgs_tibble, a list with two
tibbles will be returned. The first tibble is the orgs tibble and the second
is the uncertainty tibble which is comprised of the percent uncertainty [0-1]
of the species level assignment of each otu. If the otu was not assigned to
the species level than the uncertainty level is 1. If the otu was assigned to
only a single species than the uncertainty level is 0. If the species level
was assigned to N possible species, than the uncertainty level is 1-(1/N).

If an object of class keggerator is passed to orgs_tibble than the same two
tbls are returned, but rather than being in a list alone, they are added to
the keggerator object that is given in the orgs_tbl and species_uncert
slots, respectively.
}
\section{Methods (by class)}{
\itemize{
\item \code{tax_tbl}: method for tax_tbl

\item \code{keggerator}: method for keggerator

\item \code{phyloseq}: method for phyloseq object
}}


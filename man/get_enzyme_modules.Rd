% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-module.R
\name{get_enzyme_modules}
\alias{get_enzyme_modules}
\title{Query all modules related to each enzyme in a dataset}
\usage{
get_enzyme_modules(data, kegg_module = NULL)
}
\arguments{
\item{data}{a tibble that is output from get_pathway_enzymes()}

\item{kegg_module}{a kegg_tbl with the columns module and module_id. This
tbl can be generated with the get_kegg_module function. If NULL (default) the
KEGGerator::kegg_modules dataset will be used.}
}
\description{
Query all modules related to each enzyme in a dataset
}

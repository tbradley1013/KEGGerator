# functions that will see if a dataset has all functions required to complete
# a specified process of a certain pathway

#' returns genomes that have all the functions in specified module
#'
#' @param data a tibble that is output from \code{get_genome_orthologies()}
#' @param function_data a tibble containing a function id for each of
#' the orthologies in a desired pathway. An example of this can be seen in the
#' \code{nitrogen_metabolism_orth_fun} dataset
#' @param complete_functions a tibble containing all of the functions
#' required for each full process in a given pathway. An example can be seen
#' in the \code{nitrogen_metabolism_fun_groups} dataset
#'
#' @details This function takes three tibbles. The first (\code{data}) is a tibble that
#' is output from \code{get_genome_orthologies}. This contains all of the
#' orthologies associated with both the genome and the desired pathway.
#'
#' The second tibble (\code{function_data}) contains the function id and function meaning for each
#' orthology in a given pathway. This data should be stored in columns named
#' \code{func_id} and \code{func_meaning}, respectively. There should also be a
#' column containing the modules that these functions are associated with by
#' the name of \code{modules_name}.
#'
#' The third thibble (\code{complete_functions}) should be a tibble with two columns.
#' The first column should be named \code{modules_name} and there should be
#' a row for every module in a specific pathway. The second column should be
#' named \code{funcs_in_module} and should be a list column where each row contains
#' a list with all of the function ids (corresponding to the \code{func_id}
#' column in the function_data dataset)
#'
#' @export
has_complete_function <- function(data, function_data, complete_functions){
  if (!"tbl_df" %in% class(data)) stop("data must be a tibble output from get_genomes_orthologies()")
  if ((!"orthology_id" %in% colnames(data)) | (!"orthology_id" %in% colnames(function_data))){
    stop("both data and function_data must have a column named 'orthology'")
  }

  output <- data %>%
    dplyr::left_join(function_data, by = "orthology_id") %>%
    dplyr::group_by(genome, modules_name) %>%
    tidyr::nest() %>%
    dplyr::mutate(genome_funs = purrr::map(data, ~.x$func_id)) %>%
    dplyr::right_join(complete_functions, by = "modules_name") %>%
    dplyr::mutate(
      has_all = purrr::map2_lgl(
        genome_funs,
        funcs_in_module,
        ~{
          all((.y %in% .x))
        }
      )
    ) %>%
    dplyr::filter(has_all) %>%
    tidyr::unnest(data) %>%
    dplyr::select(-has_all)

  return(output)
}

#' find all genomes that complete all of the desired functions
#'
#' @param data a tibble containing all of the orthologies for each genome
#' in a specific pathway. This should be the output of the
#' \code{get_genome_orthologies} function
#' @param function_data a tibble containing a function id for each of
#' the orthologies in a desired pathway. An example of this can be seen in the
#' \code{nitrogen_metabolism_orth_fun} dataset
#' @param desired_functions a list of vector containing all of the function ids
#' (corresponding to those in the \code{func_id} column of the function_data
#' tibble) that a genome should carry out
#'
#' #' @details This function takes two tibbles and a vector. The first (\code{data}) is a tibble that
#' is output from \code{get_genome_orthologies}. This contains all of the
#' orthologies associated with both the genome and the desired pathway.
#'
#' The second tibble (\code{function_data}) contains the function id and function meaning for each
#' orthology in a given pathway. This data should be stored in columns named
#' \code{func_id} and \code{func_meaning}, respectively. There should also be a
#' column containing the modules that these functions are associated with by
#' the name of \code{modules_name}.
#'
#' The third argument (\code{desired_functions}) takes a vector or a list of
#' function ids (corresponding to thos in the \code{func_id} column in the
#' \code{function_data} tibble).
#'
#' @export
has_specific_function <- function(data, function_data, desired_functions){

  output <- data %>%
    dplyr::left_join(function_data, by = "orthology_id") %>%
    dplyr::filter(!is.na(func_meaning)) %>%
    dplyr::group_by(genome, modules_name) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      genome_funs = purrr::map(data, ~.x$func_id),
      has_all = purrr::map_lgl(genome_funs, ~{
        all(desired_functions %in% .x)
      })
    ) %>%
    dplyr::filter(has_all) %>%
    dplyr::select(-has_all) %>%
    tidyr::unnest(data)

  return(output)
}

#' Genomic Functional Heatmap
#'
#' This function will crteate a heatmap, either static (via ggplot2) or
#' interactive (via plotly), that will allow for the visualization of
#' genomes in a given dataset that complete different functions
#'
#'
genomic_func_heatmap <- function(data, otu_tibble = NULL, tax_tibble = NULL,
                              ps_orig = NULL, interactive = TRUE,
                              relative_abundance = FALSE,
                              abundance_threshold = NULL,
                              relative_abundance_threshold = NULL){
  if (!all(c("func_id", "func_meaning") %in% colnames(data))){
    stop("`data` must be the output of either `has_complete_function`",
         " or `has_specific_function`", call. = FALSE)
  }

  if (is.null(otu_tibble) & is.null(tax_tibble) & is.null(ps_orig)) {
    stop("either `otu_tibble` and `tax_tibble` or `ps_orig` must be supplied", call. = FALSE)
  }

  # if (is.null(otu_tibble) + is.null(ps_orig) != 1){
  #   stop("only one of `otu_tibble` and `ps_orig` may be specified")
  # }

  if (is.null(ps_orig)){
    if (is.null(otu_tibble) + is.null(tax_tibble) > 0){
      stop("if `ps_orig` is not supplied than both `otu_tibble`",
           " and `tax_tibble` must be supplied")
    }
  }

  if (is.null(otu_tibble)) otu_tibble <- otu_tibble(ps_orig)
  if (is.null(tax_tibble)) tax_tibble <- tax_tibble(ps_orig)


  otu_cols <- otu_tibble %>%
    dplyr::select(-seq) %>%
    colnames() %>%
    purrr::map(dplyr::quo_name)


  otu_tibble_long <- otu_tibble %>%
    tidyr::gather(key = sample, value = count, !!!otu_cols) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(total_count = sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(perc_count = (count/total_count)*100)  %>%
    dplyr::left_join(
      genomes_tibble(tax_tibble, drop_taxa = FALSE),
      by = "seq"
    ) %>%
    dplyr::filter(!is.na(genome))


  genomes_func_otu <- data %>%
    dplyr::left_join(otu_tibble_long, by = "genome") %>%
    dplyr::group_by(genome, modules_name) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      orthology = purrr::map_chr(data, ~paste(.x$orthology_name %>% unique(), collapse = "\n")),
      species = purrr::map_chr(data, ~paste(.x$genome_query %>% unique() %>% head(11), collapse = "\n")),
      species_count = purrr::map_int(data, ~.x$genome_query %>% unique() %>% length()),
      genome_function = purrr::map_chr(data, ~paste(.x$func_meaning %>% unique(), collapse = "\n")),
      # function_filter = purrr::map(data, ~.x$func_meaning %>% unique() %>% as.list()),
      samples = purrr::map(data, ~.x %>% dplyr::select(sample, count, perc_count))
    ) %>%
    tidyr::unnest(samples, .drop = FALSE) %>%
    dplyr::select(-data) %>%
    dplyr::filter(!is.na(species)) %>%
    dplyr::distinct()

  if (relative_abundance) {
    plot_data <- genomes_func_otu %>%
      dplyr::mutate(z = perc_count)

    legend_title <- "Relative Abundance (%)"
  } else {
    plot_data <- genomes_func_otu %>%
      dplyr::mutate(z = count)

    legend_title <- "Abundance"
  }

  if (!is.null(abundance_threshold)) {
    if (!is.numeric(abundance_threshold)) stop("`abundance_threshold` must be of class numeric")
    plot_data <- plot_data %>%
      dplyr::filter(count >= abundance_threshold)

    if (!is.null(relative_abundance_threshold)) {
      warning("Both `abundance_threshold` and `relative_abundance_threshold` were specified",
              ". Only `abundance_threshold` will be used.")
    }
  } else if (!is.null(relative_abundance_threshold)) {
    if (!is.numeric(relative_abundance_threshold)) stop("`relative_abundance_threshold` must be of class numeric")

    plot_data <- plot_data %>%
      dplyr::filter(perc_count >= relative_abundance_threshold)
  }

  # browser()


  if (interactive) {
    num_genomes <- plot_data %>%
      dplyr::distinct(genome) %>%
      nrow()

    output <- plot_data %>%
      plotly::plot_ly(
        x = ~genome,
        y = ~sample,
        z = ~z,
        type = "heatmap",
        text = ~paste(
          "<b>Genome:</b><i>", genome,
          "</i><br><b>Species:</b><i>", species,
          "</i><br><b>Species Count:</b><i>", species_count,
          "</i><br><b>Sample:</b><i>", sample,
          "</i><br><b>Abundance:</b><i>", count,
          "</i><br><b>Relative Abundance (%):</b><i>", perc_count,
          "</i><br><b>Orthologies:</b><i>", orthology,
          "</i><br><b>Functionalities:</b><i>", genome_function,
          "</i><br><b>Module:</b><i>", modules_name, "</i>"
        ),
        hoverinfo = "text",
        hoverlabel = list(
          bgcolor = "white",
          font = list(
            color = "black"
          ),
          bordercolor = "black"
        ),
        colors = grDevices::colorRamp(c("blue", "red")),
        colorbar = list(
          title = legend_title
        )
      ) %>%
      plotly::layout(
        margin = list(
          b = if_else(num_genomes > 10, 175, 100),
          l = 80
        ),
        font = list(
          size = 10
        ),
        xaxis = list(
          title = ""
        )
      )
  } else {
    output <- plot_data %>%
      ggplot2::ggplot(aes(x = genome, y = sample, fill = z)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(high = "red", low = "blue", name = legend_title) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  }


  return(output)

}

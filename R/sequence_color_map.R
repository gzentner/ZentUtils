
#' Color Map
#'
#' @description
#' Generates a tile plot of sequences of interest wherein each base is represented
#' by a specific color.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual theme theme_minimal element_blank
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#'
#' @param zu_obj A ZentUtils object
#' @param cols Colors to use for base representation. Use NA for default colors
#'             or provide a character vector of length 4.
#'
#' @return A sequence color map
#' @export
#'
#' @examples
#' color_map(zent, cols = c("green", "blue", "yellow", "red"))

color_map <- function(zu_obj, cols = NA) {

  # Input check
  if (!is(zu_obj, "zu_obj"))
    stop("input must be a ZentUtils object")

  if (length(cols) != 4 && !is.na(cols)) {
    stop("cols must be NA or a character vector of length 4")
  }

  # Generate long table of bases for plotting
  seq_length <- unique(stringr::str_length(zu_obj@seqs))
  seqs_split <- zu_obj@seqs %>%
    as.character %>%
    stringr::str_split(pattern = "", simplify = T) %>%
    as.data.frame %>%
    magrittr::set_names(1:seq_length) %>%
    tibble::rowid_to_column(var = "sequence") %>%
    tidyr::pivot_longer(-sequence, names_to = "position", values_to = "base") %>%
    dplyr::mutate(position = as.integer(position))

  ifelse(
    is.na(cols),
    base_cols <- c(
      "A" = "#109649",
      "C" = "#255C99",
      "G" = "#F7B32C",
      "T" = "#D62839"
    ),
    base_cols <- c(
      "A" = cols[1],
      "C" = cols[2],
      "G" = cols[3],
      "T" = cols[4]
    )
  )

  # Generate color plot
  p <- ggplot(seqs_split, aes(x = .data$position, y = .data$sequence)) +
    geom_tile(aes(fill = .data$base)) +
    theme_minimal() +
    scale_fill_manual(values = base_cols) +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank()
    )

  return(p)

}

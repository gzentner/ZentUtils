
#' Expand peak summits to a desired width
#'
#' @param zu_obj A ZentUtils object
#' @param length Number of bases to expand the summit
#'
#' @return A GRanges object in slot 'regions' (overwrites original regions)
#' @export
#'
#' @examples expand_summits(zu_obj, length = 100)

expand_summits <- function(zu_obj, length = 100) {

  # Check chromosome style - must be UCSC (with the "chr" prefix) to be
  # compatible with BSgenome.
  if (GenomeInfoDb::seqlevelsStyle(GenomeInfoDb::seqlevels(zu_obj@regions))[1] != "UCSC") {
      GenomeInfoDb::seqlevelsStyle(GenomeInfoDb::seqlevels(zu_obj@regions)) <- "UCSC" }

  # Expand regions. Note that 1 is subtracted from 'length' to account for
  # 0-based coordinates
  zu_obj@expanded_regions <- plyranges::stretch(zu_obj@regions, length - 1)

  return(zu_obj)
}

####

#' Retrieve sequences corresponding to ranges of interest
#'
#' @param zu_obj A ZentTools object
#' @param region_type Whether to extract sequences from original regions
#'                    ("imported") or expanded regions ("exported")
#' @param genome A BSgenome object for the organism of interest
#'
#' @return A DNAStringSet object in slot 'seqs'
#' @export
#'
#' @examples
#' saccer3 <- BSgenome.Scerevisiae.UCSC.sacCer3
#' get_seq(zu_obj, genome = saccer3)

get_seqs <- function(zu_obj,
                     region_type = "imported",
                     genome) {

  if(region_type == "imported") {seq_regions <- zu_obj@regions}
  if(region_type == "expanded") {seq_regions <- zu_obj@expanded_regions}

  zu_obj@seqs <- BSgenome::getSeq(genome, seq_regions)

  names(zu_obj@seqs) <- seq_regions$name

  return(zu_obj)

  }

#' Generate sequence color map
#'
#' @param zu_obj A ZentUtils object
#' @param cols Colors to use for base representation. Use NA for default colors
#'             or provide a character vector of length 4.
#'
#' @return A sequence color map
#' @export
#'
#' @examples color_map(zu_obj, cols = c("green", "blue", "yellow", "red"))

color_map <- function(zu_obj,
                      cols = NA) {
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
    setNames(1:seq_length) %>%
    tibble::rowid_to_column(var = "sequence") %>%
    tidyr::pivot_longer(-sequence, names_to = "position", values_to = "base")

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
  ggplot2::ggplot(seqs_split, aes(x = position, y = sequence)) +
    geom_tile(aes(fill = base)) +
    theme_minimal() +
    scale_fill_manual(values = base_cols) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )

}

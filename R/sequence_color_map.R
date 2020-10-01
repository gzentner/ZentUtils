
#' Expand regions to a desired width
#'
#' @importFrom GenomeInfoDb seqlengths
#'
#' @param zu_obj A ZentUtils object
#' @param genome A BSgenome object for the organism of interest
#' @param length Number of bases to expand the region
#'
#' @return A GRanges object in slot 'expanded_regions'
#' @export
#'
#' @examples
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3
#' zent <- expand_regions(zent, genome, 100)

expand_regions <- function(zu_obj, genome, length = 100) {

  # Expand regions. Note that 1 is subtracted from 'length' to account for
  # 0-based coordinates

  zu_obj@expanded_regions <- plyranges::stretch(zu_obj@regions, length - 1)

  # Ensure that seqnames are in the proper order
  zu_obj@expanded_regions <- GenomeInfoDb::sortSeqlevels(zu_obj@expanded_regions)

  # Set seqlengths to enable trimming
  seqlengths_df <- as.data.frame(GenomeInfoDb::seqlengths(genome)) %>%
    tibble::rownames_to_column("seqname") %>%
    dplyr::filter(seqname != "chrM")

  suppressWarnings(GenomeInfoDb::seqlengths(zu_obj@expanded_regions) <- seqlengths_df[,2])

  # Trim out-of-bounds regions
  zu_obj@expanded_regions <- GenomicRanges::trim(zu_obj@expanded_regions)

  # Drop sequences under the maximum width due to trimming
  zu_obj@expanded_regions <- zu_obj@expanded_regions[BiocGenerics::width(zu_obj@expanded_regions) ==
                                                     max(BiocGenerics::width(zu_obj@expanded_regions))]

  dropped_seq_n <- length(zu_obj@regions) - length(zu_obj@expanded_regions)

  print(paste(dropped_seq_n, "sequences were removed due to being out-of-bounds."), quote = F)

  return(zu_obj)
}

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
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3
#' zent <- get_seqs(zent, region_type = "expanded", genome = genome)

get_seqs <- function(zu_obj,
                     region_type = "expanded",
                     genome) {

  if(region_type == "imported") {seq_regions <- zu_obj@regions}
  if(region_type == "expanded") {seq_regions <- zu_obj@expanded_regions}

  zu_obj@seqs <- BSgenome::getSeq(genome, seq_regions)

  #names(zu_obj@seqs) <- seq(1:length(zent@seqs))

  return(zu_obj)

  }

#' Generate sequence color map
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
    tidyr::pivot_longer(-sequence, names_to = "position", values_to = "base") %>%
    dplyr::mutate_at(vars(position), ~ as.integer(.))

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
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )

  return(p)

}

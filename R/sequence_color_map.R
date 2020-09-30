

#' Expand regions to a desired width
#'
#' @param zu_obj A ZentUtils object
#' @param sizes Path to a two-column file of chromosome sizes from UCSC (can be a URL)
#' @param both Distance to add to the start and end of each region
#' @param left Distance to add to the start of the region
#' @param right Distance to add to the end of the region
#' @param strand Define left and right based on strand
#' @param trim Adjust coordinates for out-of-bounds intervals
#'
#' @return A GRanges object in slot 'expanded_regions'
#' @export
#'
#' @examples reb1 <- expand_regions(reb1,
#' sizes = BSgenome.Scerevisiae.UCSC.sacCer3, both = 50)

expand_regions <- function(zu_obj,
                           sizes,
                           both = 0,
                           left = 0,
                           right = 0,
                           strand = TRUE,
                           trim = TRUE) {
  # Input checks
  if (!is(zu_obj, "zu_obj"))
    stop("input must be a ZentUtils object")

  # Make both and left/right mutually exclusive
  if (both != 0 & left != 0 | right != 0)
    stop("Specify a value only for 'both' or 'left' and/or 'right'")

  # Check chromosome style - must be UCSC (with the "chr" prefix) to be
  # compatible with valr and BSgenome.
  if (GenomeInfoDb::seqlevelsStyle(GenomeInfoDb::seqlevels(zu_obj@regions))[1] != "UCSC") {
      GenomeInfoDb::seqlevelsStyle(GenomeInfoDb::seqlevels(zu_obj@regions)) <- "UCSC"
  }

  # Convert GRanges object to a BED-format dataframe for valr
  regions_df <- data.frame(
    chrom = GenomicRanges::seqnames(zu_obj@regions),
    start = BiocGenerics::start(zu_obj@regions) - 1,
    end = BiocGenerics::end(zu_obj@regions),
    name = zu_obj@regions$name,
    score = c(rep(".", length(zu_obj@regions))),
    strand = BiocGenerics::strand(zu_obj@regions)
  )

  # Read chromosome sizes file for valr
  sizes <- GenomeInfoDb::seqlengths(sizes) %>%
    as.data.frame %>%
    tibble::rownames_to_column("chrom") %>%
    rename("." = "size") %>%
    as_tibble

  # Expand regions
  if (both != 0) {
    expanded_regions <-
      valr::bed_slop(
        regions_df,
        sizes,
        both = both,
        strand = strand,
        trim = trim
      )

    zu_obj@expanded_regions <- GenomicRanges::makeGRangesFromDataFrame(expanded_regions)

  } else if (left > 0 | right > 0) {
    expanded_regions <-
      valr::bed_slop(
        regions_df,
        sizes,
        left = left,
        right = right,
        strand = strand,
        trim = trim
      )

    zu_obj@expanded_regions <- GenomicRanges::makeGRangesFromDataFrame(expanded_regions)

  }

  return(zu_obj)
}

####

#' Retrieve sequences corresponding to ranges of interest
#'
#' @param zu_obj A ZentTools object
#' @param region_type Whether to extract sequences from original regions
#'                    ("imported") or expanded regions ("expanded")
#' @param genome A BSgenome object for the organism of interest
#'
#' @return A DNAStringSet object in slot 'seqs'
#' @export
#'
#' @examples
#' saccer3 <- BSgenome.Scerevisiae.UCSC.sacCer3
#' reb1 <- get_seq(reb1, region_type = "expanded", genome = saccer3)

get_seqs <- function(zu_obj,
                     region_type = "expanded",
                     genome) {

  # Input checks
  if (!is(zu_obj, "zu_obj"))
    stop("input must be a ZentUtils object")

  if(region_type != "imported" | "expanded")
    stop("region_type must be 'imported' or 'expanded'")

  if (region_type == "imported") {
    seq_regions <- zu_obj@regions
  }
  if (region_type == "expanded") {
    seq_regions <- zu_obj@expanded_regions
  }

  zu_obj@seqs <- BSgenome::getSeq(genome, seq_regions)

  names(zu_obj@seqs) <- seq_regions$name

  return(zu_obj)

}

#' Generate sequence color map
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual theme element_blank
#' @importFrom magrittr %>%
#'
#' @param zu_obj A ZentUtils object
#' @param cols Colors to use for base representation. Use NA for default colors
#'             or provide a character vector of length 4.
#'
#' @return A sequence color map
#' @export
#'
#' @examples
#' color_map(reb1, cols = c("green", "blue", "yellow", "red"))

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
    stats::setNames(1:seq_length) %>%
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
  p <-
    ggplot(seqs_split, aes(x = .data$position, y = .data$sequence)) +
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


#' Sort Regions
#'
#' @description
#' Sorts imported regions by chromosomal location or score.
#'
#' @importFrom plyranges anchor_center
#'
#' @param zu_obj A ZentUtils object.
#' @param by Sort by chromosome + coordinates + start ("coord") or column 5
#'           of the imported ranges (automatically named "score" upon import).
#' @param decreasing Whether the GRanges object should be sorted in descending order.
#'
#' @return A sorted GRanges object in the 'regions' slot of the ZentUtils object
#'         (overwrites original region order).
#' @export
#'
#' @examples
#' zent <- sort_regions(zent, by = "coord", decreasing = FALSE)

sort_regions <- function(zu_obj, by = "coord", decreasing = FALSE) {

  # Input checks
  if (!is(zu_obj, "zu_obj"))
    stop("Input must be a ZentUtils object.")

  if(by == "coord") {
    zu_obj@regions <- BiocGenerics::sort(zu_obj@regions, decreasing = decreasing)
  } else if (by == "score") {
    zu_obj@regions <- BiocGenerics::sort(zu_obj@regions, by = ~score, decreasing = decreasing)
  } else {
    stop("Regions must be sorted by 'coord' or 'score'.")
  }

  return(zu_obj)

}

#' Expand Regions
#'
#' @description
#' Expands regions to a desired width. Sequences that are partially out-of-bounds
#' after expansion are removed.
#'
#' @param zu_obj A ZentUtils object.
#' @param length Number of bp to add to both sides of each region.
#'
#'
#' @return A GRanges object in the 'expanded_regions' slot of the ZentUtils object.
#' @export
#'
#' @examples
#' zent <- expand_regions(zent, length = 50)

expand_regions <- function(zu_obj, length = 10) {

  # Expand regions
  zu_obj@expanded_regions <- suppressWarnings(
    zu_obj@regions %>%
    plyranges::anchor_center() %>%
    plyranges::stretch(length * 2)
  )

  # Trim out-of-bounds regions
  zu_obj@expanded_regions <- GenomicRanges::trim(zu_obj@expanded_regions)

  # Drop sequences under the maximum width due to trimming
  zu_obj@expanded_regions <- zu_obj@expanded_regions[BiocGenerics::width(zu_obj@expanded_regions) ==
                                                     max(BiocGenerics::width(zu_obj@expanded_regions))]

  dropped_seq_n <- length(zu_obj@regions) - length(zu_obj@expanded_regions)

  message(paste(dropped_seq_n, "sequences were removed due to being out-of-bounds."))

  # Store extension length as metadata for color map labeling
  metadata(zu_obj@expanded_regions)$length <- length

  return(zu_obj)

}

#' Get sequences
#'
#' @description
#' Retrieve sequences corresponding to ranges of interest
#'
#' @param zu_obj A ZentTools object.
#' @param region_type Whether to extract sequences from original regions
#'                    ("imported") or expanded regions ("exported").
#' @param genome A BSgenome object for the organism of interest.
#'
#' @return A DNAStringSet object in the 'seqs' slot of the ZentUtils object.
#' @export
#'
#' @examples
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3
#' zent <- get_seqs(zent, region_type = "expanded", genome = genome)

get_seqs <- function(zu_obj, region_type = "expanded", genome) {

  if(region_type == "imported") {seq_regions <- zu_obj@regions}
  if(region_type == "expanded") {seq_regions <- zu_obj@expanded_regions}

  zu_obj@seqs <- BSgenome::getSeq(genome, seq_regions)

  names(zu_obj@seqs) <- seq(1:length(zent@seqs))

  return(zu_obj)

}

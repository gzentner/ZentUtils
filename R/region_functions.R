
#' Sort Regions
#'
#' @description
#' Sorts imported regions by chromosomal location or score.
#'
#' @importFrom plyranges anchor_center
#'
#' @param zu_obj A ZentUtils object.
#' @param by Sort by chromosomal coordinates ("coord") or column 5 of the imported
#'          ranges (automatically named "score" upon import).
#' @param decreasing Whether the GRanges object should be sorted in descending order.
#'
#' @return A sorted GRanges object in the 'regions' slot of the ZentUtils object
#'         (overwrites original region order).
#' @export
#'
#' @examples
#' zent <- sort_regions(zent, by = "score", decreasing = TRUE)

sort_regions <- function(zu_obj, by = "coord", decreasing = FALSE) {

  # Input checks
  if (!is(zu_obj, "zu_obj"))
    stop("Input must be a ZentUtils object.")

  by <- match.arg(stringr::str_to_lower(by), choices=c("coord", "score"))

  # Sort
  zu_obj@regions <- switch(by,
    "coord"=BiocGenerics::sort(zu_obj@regions, decreasing = decreasing),
    "score"=BiocGenerics::sort(zu_obj@regions, by = ~score, decreasing = decreasing)
  )

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
#' zent <- expand_regions(zent, length = 10)

expand_regions <- function(zu_obj, length = 10) {

  # Input checks
  if (!is(zu_obj, "zu_obj")) stop("zu_obj must be a 'zu_obj'")
  if (!is(length, "numeric") || length < 1 || length %% 1 != 0) {
    stop("length must be a positive integer greater than or equal to 1")
  }

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
#'                    ("imported") or expanded regions ("expanded").
#' @param genome A BSgenome object for the organism of interest.
#'
#' @return A DNAStringSet object in the 'seqs' slot of the ZentUtils object.
#' @export
#'
#' @examples
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3
#' zent <- get_seqs(zent, region_type = "expanded", genome = genome)

get_seqs <- function(zu_obj, region_type = "expanded", genome) {

  # Input checks.
  if (!is(zu_obj, "zu_obj")) stop("zu_obj must be a 'zu_obj'")
  region_type <- match_arg(
    stringr::str_to_lower(region_type),
    choices=c("imported", "expanded")
  )

  # Extract proper regions.
  seq_regions <- switch(region_type,
    "imported"=zu_obj@regions,
    "expanded"=zu_obj@expanded_regions
  )

  # Retrieve the sequences.
  zu_obj@seqs <- BSgenome::getSeq(genome, seq_regions)
  names(zu_obj@seqs) <- seq_len(length(zent@seqs))

  return(zu_obj)

}

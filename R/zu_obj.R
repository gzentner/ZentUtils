
#' ZentUtils class definition
#'
#' @slot regions A GRanges object containing imported regions
#' @slot expanded_regions A GRanges object containing expanded regions
#' @slot seqs A DNAStringSet object containing sequences of regions of interest
#'
#' @rdname ZentUtils-class
#'
#' @export

setClass(
  "zu_obj",
  slots = c(
    regions = "GRanges",
    expanded_regions = "GRanges",
    seqs = "DNAStringSet"
    )
)


#' ZentUtils constructor function
#'
#' @import methods
#'
#' @param data A BED file containing the regions of interest
#'
#' @return A ZentUtils object
#' @export
#'
#' @examples
#' reb1 <- zentutils(system.file("extdata", "reb1_motifs.bed", package = "ZentUtils"))

zentutils <- function(data) {

  zu_obj <- new(
    "zu_obj",
    regions = rtracklayer::import(data)
  )

  return(zu_obj)
}


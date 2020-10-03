
#' ZentUtils class definition
#'
#' @slot regions A GRanges object containing imported regions.
#' @slot expanded_regions A GRanges object containing expanded regions.
#' @slot seqs A DNAStringSet object containing sequences of regions of interest.
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

#' ZentUtils Constructor Function
#'
#' @description
#' Creates a new ZentUtils object that stores a GRanges object generated from
#' an imported BED file.
#'
#' @import methods
#'
#' @param data A BED file containing the regions of interest.
#' @param genome Identifier for the genome from which the regions were derived.
#'
#' @details
#' \itemize{
#' \item{The BED file should contain 6 columns, with chromosome, start, and
#' end coordinates in columns 1-3 and strand in column 6.}
#' \item{BED columns names are converted to "chrom", "start", "end",
#' "region_name", "score", and "strand" upon import.}
#' }
#'
#' @return A new ZentUtils object with a GRanges object in the 'regions' slot.
#' @export
#'
#' @examples
#' zent <- zentutils(system.file("extdata", "homer_reb1_badis_motif_scan.bed", package = "ZentUtils"),
#'                   genome = "sacCer3")

zentutils <- function(data, genome = NA) {

  zu_obj <- new(
    "zu_obj",
    regions = readr::read_tsv(data, col_names = c("chrom", "start", "end",
                                                  "region_name", "score", "strand")) %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  )

  # If necessary, add "chr" prefix to ENSEMBL or NCBI chromosome names
  if (GenomeInfoDb::seqlevelsStyle(GenomeInfoDb::seqlevels(zu_obj@regions))[1] != "UCSC") {
      GenomeInfoDb::seqlevelsStyle(GenomeInfoDb::seqlevels(zu_obj@regions)) <- "UCSC" }

  # Remove any mitochondrial regions
  zu_obj@regions <- GenomeInfoDb::dropSeqlevels(zu_obj@regions, value = "chrM", pruning.mode = "coarse")

  # Load list of chromosome length vectors
  genome_list <- readRDS("data/seqlengths.rds")

  # Add chromosome lengths and genome to GRanges object in the 'regions' slot of the new ZentUtils object
  GenomeInfoDb::seqlengths(zu_obj@regions) <- as.vector(genome_list[[genome]])
  GenomeInfoDb::genome(zu_obj@regions) <- genome

  return(zu_obj)
}


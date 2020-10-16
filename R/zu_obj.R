
#' ZentUtils class definition
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings DNAStringSet
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
  slots = list(
    regions = "GRanges",
    expanded_regions = "GRanges",
    seqs = "DNAStringSet"
    ),
  prototype = list(
    regions = GRanges(),
    expanded_regions = GRanges(),
    seqs = DNAStringSet()
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
#' @param mito_chr Seqname of mitochondrial chromosome.
#'
#' @details
#' \itemize{
#' \item{The BED file should contain 6 columns, with chromosome, start, and
#' end coordinates in columns 1-3 and strand in column 6.}
#' \item{BED columns names are converted to "seqnames", "start", "end",
#' "region_name", "score", and "strand" upon import.}
#' }
#'
#' @return A new ZentUtils object with a GRanges object in the 'regions' slot.
#' @export
#'
#' @examples
#' zent <- zentutils(system.file("extdata", "homer_reb1_badis_motif_scan.bed", package = "ZentUtils"), genome = "sacCer3", mito_chr = "chrM")

zentutils <- function(data, genome = NA, mito_chr = "chrM") {

  # Read in data.
  bed <- data %>%
    readr::read_tsv(col_names = c("seqnames", "start", "end", "region_name", "score", "strand")) %>%
    plyranges::as_granges()

  # Create a new zu_obj.
  zu_obj <- new(
    "zu_obj",
    regions = bed
  )

  # If necessary, add "chr" prefix to ENSEMBL or NCBI chromosome names
  if (GenomeInfoDb::seqlevelsStyle(GenomeInfoDb::seqlevels(zu_obj@regions))[1] != "UCSC") {
      GenomeInfoDb::seqlevelsStyle(GenomeInfoDb::seqlevels(zu_obj@regions)) <- "UCSC" }

  # Sort seqlevels to ensure compatibility with seqinfo
  zu_obj@regions <- GenomeInfoDb::sortSeqlevels(zu_obj@regions)

  # Retrieve genome information from UCSC
  genome_info <- plyranges::genome_info(genome = genome)

  zu_obj@regions <- zu_obj@regions %>%
    plyranges::set_genome_info(genome = genome,
                               seqnames = genome_info@seqinfo@seqnames,
                               seqlengths = genome_info@seqinfo@seqlengths,
                               is_circular = genome_info@seqinfo@is_circular)

  # Remove any mitochondrial regions
  zu_obj@regions <- GenomeInfoDb::dropSeqlevels(zu_obj@regions, value = mito_chr, pruning.mode = "coarse")

  return(zu_obj)
}


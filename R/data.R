#' Chromosome lengths
#'
#' @description
#' Chromosome names and lengths for commonly used reference genomes.
#'
#' @format A list of named numeric vectors for the following genome assemblies:
#' \itemize{
#'   \item{hg19}
#'   \item{hg38}
#'   \item{mm9}
#'   \item{mm10}
#'   \item{dm6}
#'   \item{sacCer3}
#' }
#' @source
#' Chromosome names and lengths were extracted from BSgenome objects using the following code template:
#'
#' \code{seqlengths[["sacCer3"]] <- setNames(as.vector(
#' seqlengths(BSgenome.Scerevisiae.UCSC.sacCer3)[1:16]),
#' as.vector(seqnames(BSgenome.Scerevisiae.UCSC.sacCer3)))}
"seqlengths"

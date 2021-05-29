#' Compute sequence correlation matrix with kernel K_hat^3
#'
#' @description Uses computing kernel to compute sequence similarity. Similarity
#' is expressed as a number from range [0, 1], where 1 means equality.
#'
#' @param x [\code{character}]\cr
#'  A vector of sequences this function is applied to.
#' @param max_kmer_length [\code{integer(1)}]\cr
#'  Maximum k-mer length used for computations. Higher value results in longer
#'  computation time, but gives more accurate scores. Defaults to 5.
#' @param exponential [\code{numeric(1)}]\cr
#'  A parameter used to control impact of differences between sequences. Higher
#'  value usually results in lower scores. Default value is 0.125, as suggested
#'  in the paper.
#'
#' @return A symmetric numeric matrix with dimensions n x n, where n is the number of
#' sequences in supplied sequence vector. A value under index (i, j) is equal to the
#' similarity of sequences i and j.
#'
#' @details
#' This algorithm is based on the paper called "Towards a Mathematical
#' Foundation of Immunology and Amino Acid Chains" by Wen-Jun Shen,
#' Hau-San Wong, Quan-Wu Xiao, Xin Guo and Stephen Smale (arXiv:1205.6031).
#'
#' This sequence correlation method does not employ any alignment methods in its
#' computations. Instead, it bases heavily on k-mers, comparing analogous k-mers
#' across sequences.
#'
#' @examples
#' sq_ami <- c("PPAVMMFDILKKIQ", "PQEWYTWLPVMCTN", "PQKWLANMMAQ")
#' result <- computing_kernel(sq_ami)
#'
#' # Results in a symmetric matrix 3x3:
#' result
#'
#' # Accessing similarity between sequences number 3 and 1:
#' result[3, 1]
#' result[1, 3]  # the same value as above
#'
#' # Using higher exponential results amplifies differences:
#' computing_kernel(sq_ami, exponential = 0.5)
#'
#' # max_kmer_length can be higher than sequence length.
#' # In that case longer k-mers are simply omitted:
#' computing_kernel(sq_ami, max_kmer_length = 25)
#'
#' @export
computing_kernel <- function(x, max_kmer_length = 5, exponential = 0.125, ...)
  UseMethod("computing_kernel")

#' @export
computing_kernel.default <- function(x, max_kmer_length = 5, exponential = 0.125, ...)
  stop("'computing_kernel' isn't implemented for this type of object", call. = FALSE)

#' @rdname computing_kernel
#' @export
computing_kernel.character <- function(x, max_kmer_length = 5, exponential = 0.125, ...) {
  assert_count(max_kmer_length)
  assert_number(exponential, lower = 0)
  
  CPP_computing_kernel(x, max_kmer_length, exponential)
}

#' MutSim: A package for generating cancer mutations in genes according to mutational processes.
#'
#' The MutSim package allows to generate mutations in a set of genes according to a
#' given mutational process, total number of silent mutations per gene, total numbers
#' of silent mutations per sample and total number of non-silent mutations per sample.
#' It allows to observe the distribution of non-silent mutations in the genes of
#' interest and estimate the probability of damaging mutations occuring by chance.
#' Main function, \code{generate_mutation}, takes a set of genes, list of samples,
#' respective numbers of silent/nonsilent mutations observed, and a mutation distribution
#' of a mutagenic process, and generates a set of mutations if the genes of interest.
#'
#' @section MutSim functions:
#' generate_mutations()
#' annotate_simulated_variants()
#'
#' @docType package
#' @name MutSim
NULL

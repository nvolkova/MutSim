#' @include gene_features_class.R

# utility functions

##### --------------------------------------------------------------------------------

#' Matching HUGO SYMBOL gene names with TxDb gene names
#'
#' @param geneset Set of genes
#' @param txdb A TxDb object
#' @return Vector of TxDb names
#'
TxDbnames <- function(geneset, txdb) {

  hg19.genes <- genes(txdb)

  gene_symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = hg19.genes$gene_id,
                                       columns = "SYMBOL",
                                       keytype = "ENTREZID")

  hg19.genes$gene_id <- gene_symbol$SYMBOL

  inds <- names(hg19.genes)[match(geneset, hg19.genes$gene_id)]

  TxDb_genes <- tryCatch(
    {
      names(hg19.genes[inds])
    },
    error=function(cond) {
      message(paste("Some gene names were not found: ", geneset[is.na(inds)]))
      message("Please ensure you are using correct gene names.")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("Name conversion caused a warning:", url))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    })

  return(TxDb_genes)
}

##### --------------------------------------------------------------------------------

#' Matching TxDb gene names with their chromosomes
#'
#' @param TxDb_genes Set of TxDb gene names
#' @param transcripts A list of GRanges of exons per gene
#' @return Vector of chromosome names

get_chromosomes <- function(TxDb_genes, transcripts) {

  chrom <- sapply(TxDb_genes, function(x) {

    return(as.character(seqnames(transcripts[[x]]))[1])

  })

  return(chrom)

}

##### --------------------------------------------------------------------------------


#' Simplifying GRanges
#'
#' @description  Aggregate intersecting exons and select the smallest non-intersecting
#'               set of ranges.
#' @param gr GRanges object with ranges of exons
#' @return Reduced GRanges object with ranges of exons
min_ranges <- function(gr) {

  if (length(gr) == 1) return(gr)

  gr_start <- start(gr)

  gr_end <- end(gr)

  to.merge <- vector('numeric',length(gr))

  m = 1

  for (j in 2:length(gr)) {

    if (gr_start[j] < gr_end[j-1]) {

      if (gr_end[j] <= gr_end[j-1]) {

        to.merge[j] <- -1

      } else {

        to.merge[j] <- m

        to.merge[j-1] <- m

        m <- m+1

      }

    }

  }

  if (sum(to.merge==-1)>0) {

    gr <- gr[-which(to.merge == -1)]

    to.merge <- to.merge[-which(to.merge == -1)]

  }

  if (sum(to.merge>0)>0) {

    for (j in 1:max(to.merge)) {

      start(gr[to.merge == j]) <- min(gr_start[to.merge == j])

      end(gr[to.merge == j]) <- min(gr_end[to.merge == j])

    }

  }

  return(unique(gr))

}

##### --------------------------------------------------------------------------------

#' Simple reverse complement for mutation types
#'
#' @param nucl Mutation type name, e.g. \code{"A[C>G]B"}
#' @return Reverse complement of it
RevComMutType <- function(nucl) {

  tmp <- unlist(strsplit(nucl, split = ''))

  return(paste0(RevCom(tmp[7]),'[',RevCom(tmp[3]),'>',RevCom(tmp[5]),']',RevCom(tmp[1])))

}
#' Simple reverse complement for trinuleotides
#'
#' @param nucl Character, trinucletide sequence, e.g. 'ATA'
#' @return Reverse complement of it
RevCom <- function(x) {

  return(as.character(reverseComplement(DNAString(x))))

}

##### --------------------------------------------------------------------------------

#' Similarity function for signatures
#'
#' @param x,y Vectors of mutation probability of same length
#' @return Cosine similarity score
#' @references [Cosine similarity](https://en.wikipedia.org/wiki/Cosine_similarity)
cosine <- function(x,y) {

  x %*% y / sqrt(sum(x**2)) / sqrt(sum(y**2))

}

##### --------------------------------------------------------------------------------


#' Calculate mutation probability distribution per gene
#'
#' @description Adjusts the signature mutation probabilities to the ratio of trinucleotide context
#'              in the gene to that across the exome, which COSMIC mutation set corresponds to.
#' @param gene_features List of GeneFeature class objects
#' @param signature Mutation probability distribution for a mutational process
#' @return A list of numeric vectors corresponding to the mutation probability distribution across each gene.
get_all_probabilities <- function(geneset,
                                  gene_features,
                                  signature) {

  new_probs_long <- list()

  for (x in geneset) {

    new_probs <- c(signature[types.full],
                   signature[types.full])/
                      tri.counts.exome[c(types, types), 1] *
                      tri.counts(gene_features[[x]])[new.types]

    new_probs <- new_probs/sum(new_probs)

    names(new_probs) <- new.types.full

    new_probs_temp <- sapply(unique(new.types), function(y) {

      return(sum(new_probs[which(new.types == y)]/
                   tri.counts(gene_features[[x]])[y]))

    })

    new_probs_long[[x]] <- new_probs_temp[trinucleotides(gene_features[[x]])]

  }

  return(new_probs_long)

}

##### --------------------------------------------------------------------------------

#' Get the required features for each gene
#'
#' @description Supplies chromosome name, list of positions and contexts and a vector of
#'              trinucleotide counts per gene.
#' @param geneset A vector of gene names
#' @param txdb A TxDb object
#' @param genome A BSgenome object
#' @return A list of GeneFeatures objects for each gene
get_all_features <- function(geneset, txdb, genome, targets) {

  TxDb_geneset <- TxDbnames(geneset, txdb)

  names(TxDb_geneset) <- geneset

  gene_exons <- cdsBy(txdb, "gene")

  gene_exons <- gene_exons[TxDb_geneset]

  gene_features <- list()

  for (x in geneset) {

    gene_features[[x]] <- new(Class = 'GeneFeatures',
                              name = x,
                              positions = 0,
                              gen_chrom = 'NA',
                              trinucleotides = 'NA',
                              tri.counts = 0)

    gen_chrom(gene_features[[x]]) <- get_chromosomes(TxDb_geneset[x], gene_exons)

    gr <- GenomicRanges::intersect(min_ranges(gene_exons[[TxDb_geneset[x]]]),targets,ignore.strand=TRUE)

    seq <- getSeq(genome, gr)

    tri.counts(gene_features[[x]]) <- rowSums(sapply(seq, function(x) trinucleotideFrequency(x)))

    trinucleotides(gene_features[[x]]) <- unlist(lapply(1:length(seq), function(j) {

      sst <- strsplit(as.character(seq[j]), "")[[1]]

      return(paste0(sst[1:(length(sst) - 2)], sst[-c(1, length(sst))], sst[-c(1:2)]))

    }))

    tmp_positions <- NULL

    for (j in 1:length(ranges(gr))) {

      tmp_positions <- c(tmp_positions, c((start(ranges(gr))[j] + 1):(end(ranges(gr))[j] - 1)))

    }

    positions(gene_features[[x]]) <- tmp_positions

    print(paste0('Features collected for ', x,', ',match(x,geneset),' out of ',length(geneset),' done.'))

  }

  return(gene_features)

}


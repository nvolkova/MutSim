
#' An S4 class to represent a set of features of a gene.
#'
#' @slot name Name of the gene
#' @slot positions Positions of its exons in the genome
#' @slot gen_chrom Chomosome name
#' @slot trinucleotides The sequence of contexts in its exons
#' @slot tri.counts Trinucleotide tri.counts per gene exons
setClass("GeneFeatures",

         slots=list(

           name="character",

           positions="numeric",

           gen_chrom="character",

           trinucleotides="character",

           tri.counts="numeric"

         ))

### Name ------------------------------------------------

#' Name generic
#'
setGeneric("name", function(object, ...) standardGeneric("name"))


#' Assignable name generic
#'
setGeneric("name<-", function(object, value) standardGeneric("name<-"))


#' Name method for GeneFeatures class
#'
#' @param x An object of class GeneFeatures
#' @return \code{name} field
#'
setMethod("name", "GeneFeatures", function(object, ...) object@name)


#' Name replacement method for GeneFeatures class
#'
#' @param x A character
#' @return An object of class GeneFeatures with new \code{name} field
#'
setReplaceMethod("name", signature(object="GeneFeatures", value="character"),
                 function(object, value){
                   object@name <- value
                   return(object)
                 })

### Chromosome ------------------------------------------------

#' Chromosome generic
#'
setGeneric("gen_chrom", function(object, ...) standardGeneric("gen_chrom"))


#' Assignable chromosome generic
#'
setGeneric("gen_chrom<-", function(object, value) standardGeneric("gen_chrom<-"))


#' Chromosome method for GeneFeatures class
#'
#' @param x A GeneFeatures class object
#' @return \code{gen_chrom}field
#'
setMethod("gen_chrom", "GeneFeatures", function(object, ...) object@gen_chrom)


#' Chromosome replacement method for GeneFeatures class
#'
#' @param x A character
#' @return An object of class GeneFeatures with new \code{gen_chrom} field
#'
setReplaceMethod("gen_chrom", signature(object="GeneFeatures", value="character"),
                 function(object, value){
                   object@gen_chrom <- value
                   return(object)
                 })

### Positions ------------------------------------------------

#' Positions generic
#'
setGeneric("positions", function(object, ...) standardGeneric("positions"))


#' Assignable positions generic
#'
setGeneric("positions<-", function(object, value) standardGeneric("positions<-"))


#' Positions method for GeneFeatures class
#'
#' @param x A numeric vector
#' @return \code{positions} field
#'
setMethod("positions", "GeneFeatures", function(object, ...) object@positions)


#' Positions replacement method for GeneFeatures class
#'
#' @param x A numeric vector
#' @return An object of class GeneFeatures with new \code{positions} field
#'
setReplaceMethod("positions", signature(object="GeneFeatures", value="numeric"),
                 function(object, value){
                   object@positions <- value
                   return(object)
                 })


### Tri.counts ------------------------------------------------


#' tri.counts generic
#'
setGeneric("tri.counts", function(object, ...) standardGeneric("tri.counts"))


#' Assignable tri.counts generic
#'
setGeneric("tri.counts<-", function(object, value) standardGeneric("tri.counts<-"))



#' tri.counts method for GeneFeatures class
#'
#' @param x A GeneFeatures class object
#' @return \code{tri.counts} field
#'
setMethod("tri.counts", "GeneFeatures", function(object, ...) object@tri.counts)


#' tri.counts replacement method for GeneFeatures class
#'
#' @param x A numeric vector
#' @return An object of class GeneFeatures with new \code{tri.counts} field
#'
setReplaceMethod("tri.counts", signature(object="GeneFeatures", value="numeric"),
                 function(object, value){
                   object@tri.counts <- value
                   return(object)
                 })


### Trinucleotides ------------------------------------------------

#' Trinucleotides generic
#'
setGeneric("trinucleotides", function(object, ...) standardGeneric("trinucleotides"))

#' Assignable trinucleotides generic
#'
setGeneric("trinucleotides<-", function(object, value) standardGeneric("trinucleotides<-"))

#' Trinucleotides method for GeneFeatures class
#'
#' @param x A GeneFeatures class object
#' @return \code{trinucleotides} field
#'
setMethod("trinucleotides", "GeneFeatures", function(object, ...) object@trinucleotides)


#' Trinucleotides replacement method for GeneFeatures class
#'
#' @param x A vector of characters
#' @return An object of class GeneFeatures with new \code{trinucleotides} field
#'
setReplaceMethod("trinucleotides", signature(object="GeneFeatures", value="character"),
                 function(object, value){
                   object@trinucleotides <- value
                   return(object)
                 })



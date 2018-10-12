#' Get BAM filenames for demo
#' 
#' @return A vector of BAM filenames
#' 
#' @examples
#' bam_files <- get_demo_bam_files()
get_demo_bam_files <- function() {
    c(system.file("extdata/align1_rep1.bam", package="metagene2"),
        system.file("extdata/align1_rep2.bam", package="metagene2"),
        system.file("extdata/align2_rep1.bam", package="metagene2"),
        system.file("extdata/align2_rep2.bam", package="metagene2"),
        system.file("extdata/ctrl.bam", package="metagene2"))
}

#' Get BAM filenames for demo
#' 
#' @return A vector of BAM filenames
#' 
#' @examples
#' bam_files <- get_demo_rna_bam_files()
get_demo_rna_bam_files <- function() {
    c(system.file("extdata/cyto4.bam", package="metagene"),
                   system.file("extdata/cyto3.bam", package="metagene2"),
                   system.file("extdata/nuc4.bam", package="metagene2"),
                   system.file("extdata/nuc3.bam", package="metagene2"))
}

#' Get regions filenames for demo
#' 
#' @return A vector of regions filenames
#' 
#' @examples
#' regions <- get_demo_regions()
get_demo_region_filenames <- function() {
    c(system.file("extdata/list1.bed", package="metagene2"),
        system.file("extdata/list2.bed", package="metagene2"))
}

#' Get demo regions
#' 
#' @return A vector of regions filenames
#' 
#' @examples
#' regions <- get_demo_regions()
get_demo_regions <- function() {
    regions_list <- lapply(get_demo_region_filenames(), rtracklayer::import, format="bed")
    regions_grl <- GRangesList(regions_list)
    names(regions_grl) <-  gsub(".bed", "", basename(get_demo_region_filenames()))
    
    # We now have a named GRangesList with two set of 50 regions.
    regions_grl    
}

#' Get demo regions
#' 
#' @return A GRangesList with two genes
#' 
#' @examples
#' regions <- get_demo_rna_regions()
get_demo_rna_regions <- function() {
    gene_files = c(system.file("extdata/DPM1.bed", package="metagene2"),
                   system.file("extdata/NDUFAB1.bed", package="metagene2"))
    regions_list <- lapply(gene_files, rtracklayer::import, format="bed")
    regions_grl <- GRangesList(regions_list)
    names(regions_grl) <- c("DPM1", "NDUFAB1")
    
    # We now have a named GRangesList with the exons for two genes.
    regions_grl    
}


#' Get a demo metagene object
#'
#' @return A metagene object
#'
#' @examples
#' mg <- get_demo_metagene()
get_demo_metagene <- function() {
    regions <- get_demo_regions()
    bam_files <- get_demo_bam_files()
    metagene2$new(regions = regions, bam_files = bam_files)
}

#' Get a demo design object
#'
#' @return A \code{data.frame} corresponding to a valid design.
#'
#' @examples
#' mg <- get_demo_design()
get_demo_design <- function() {
    return(data.frame(Samples = get_demo_bam_files(),
                      align1 = c(1,1,0,0,2),
                      align2 = c(0,0,1,1,2)))
}

get_not_indexed_bam_file <- function() {
    system.file("extdata/not_indexed.bam", package = "metagene2")
}

get_different_seqnames_bam_file <- function() {
    system.file("extdata/different_header.bam", package = "metagene2")
}

get_coverage_bam_file <- function() {
    system.file("extdata/coverage.bam", package = "metagene2")
}

get_coverage_region <- function() {
    system.file("extdata/list_coverage.bed", package = "metagene2")
}

get_narrowpeak_file <- function() {
    system.file("extdata/list1.narrowPeak", package = "metagene2")
}

get_broadpeak_file <- function() {
    system.file("extdata/list1.broadPeak", package = "metagene2")
}

get_gff_file <- function() {
    system.file("extdata/test.gff", package = "metagene2")
}

get_gtf_file <- function() {
    system.file("extdata/test.gtf", package = "metagene2")
}
# Split a coverage object by its seqnames,
# and calculate 
get_subtable = function(coverages, gr, bcount) {
    grl <- split(gr, GenomeInfoDb::seqnames(gr))
    i <- vapply(grl, length, numeric(1)) > 0
    do.call("c", lapply(grl[i], get_view_means,
                        bcount = bcount, cov = coverages))
}

get_view_means = function(gr, bcount, cov) {
    chr <- unique(as.character(GenomeInfoDb::seqnames(gr)))
    gr <- unlist(tile(unlist(gr), n=bcount))
    stopifnot(length(chr) == 1)
    views <- Views(cov[[chr]], start(gr), end(gr))
    viewMeans(views)
}

bin_contiguous_regions <- function(coverage, regions, bin_count) {
  m <-  matrix(get_subtable(coverage, regions, bin_count), ncol=bin_count, byrow=TRUE)
    
  mr <- m[,bin_count:1]
  i <-as.logical(strand(unlist(regions))=="-")
  m[i,] <- mr[i,]
  
  m
}

### Function copied and modified from GenomicFeatures::coverageByTranscript
### The original function accepted input on which coverage() could be called,
### whereas this one presumes the coverage to be already computed.
###
### Most comments and notes are left as-is.
###
### Computing coverage by transcript (or CDS) of a set of ranges is an
### operation that feels a lot like extracting transcript sequences from a
### genome. Defined as an ordinary function for now.
discontiguous_coverage <- function(cvg, transcripts)
{
  ## STEP 2 - Compute unique exons ('uex').

    ex <- unlist(transcripts, use.names=FALSE)
    ## We could simply do 'uex <- unique(ex)' here but we're going to need
    ## 'sm' and 'is_unique' later to compute the "reverse index" so we compute
    ## them now and use them to extract the unique exons. That way we hash
    ## 'ex' only once (the expensive operation).
    sm <- selfmatch(ex)  # uses a hash table internally
    is_unique <- sm == seq_along(sm)
    uex2ex <- which(is_unique)  # index of unique exons
    uex <- ex[uex2ex]  # unique exons

  ## STEP 3 - Compute coverage for each unique exon ('uex_cvg').

    uex_cvg <- cvg[uex]  # parallel to 'uex'


  ## STEP 4 - Flip coverage for exons on minus strand.

    ## It feels like this is not as fast as it could be (the bottleneck being
    ## subsetting an Rle object which needs to be revisited at some point).
    uex_cvg <- revElements(uex_cvg, strand(uex) == "-")

  ## STEP 5 - Compute coverage by original exon ('ex_cvg').

    ex2uex <- (seq_along(sm) - cumsum(!is_unique))[sm]  # reverse index
    #stopifnot(identical(ex2uex[uex2ex], seq_along(uex2ex)))  # sanity
    #stopifnot(identical(ex2uex[sm], ex2uex))  # sanity
    #stopifnot(all(uex[ex2uex] == ex))  # sanity

    ex_cvg <- uex_cvg[ex2uex]  # parallel go 'ex'

  ## STEP 6 - Compute coverage of each transcript by concatenating coverage of
  ##          its exons.

    ans <- IRanges:::regroupBySupergroup(ex_cvg, transcripts)

    ans
}

bin_rle_list = function(x, bin_count) { 
    viewMeans(Views(x, breakInChunks(length(x), nchunk=bin_count)))
}

bin_discontiguous_regions <- function(coverage, regions, bin_count) {
    rle_list = discontiguous_coverage(coverage, regions)
    bin_list = lapply(rle_list, bin_rle_list, bin_count=100)
    matrix(unlist(bin_list), ncol=bin_count, byrow=TRUE)
}

bin_region_coverages = function(coverages, regions, bin_count) {
    results = list()
    for(cov_name in names(coverages)) {
        if(is(regions, "GRangesList")) {
            results[[cov_name]] = bin_discontiguous_regions(coverages[[cov_name]], regions, bin_count)
        } else {
            results[[cov_name]] = bin_contiguous_regions(coverages[[cov_name]], regions, bin_count)
        }        
    }
    return(results)
}

bin_coverages_s = function(coverage_s, regions, bin_count) {
    if(is.null(coverage_s)) {
        return(NULL)
    } else {
        return(bin_region_coverages(coverage_s, regions, bin_count))
    }
}

group_coverages_s = function(coverage_s, design, noise_removal, bam_handler=NULL) {
    if(is.null(coverage_s)) {
        results = NULL
    } else {
        if (!is.null(noise_removal)) {
            results <- remove_controls(coverage_s, design, bam_handler)
        } else {
            results <- merge_chip(coverage_s, design)
        }
    }
    
    return(results)
}

remove_controls = function(coverages, design, bam_handler) {
    results <- list()
    for (design_name in colnames(design)[-1]) {
        # Add up coverage for all ChIP and all input bams.
        chip_results <- merge_reduce(coverages, design, design_name, 1)
        input_results <- merge_reduce(coverages, design, design_name, 2)
        
        
        if (length(input_results$BamNames) > 0) {
            # If we had input bams, perform noise reduction.
            noise_ratio <-
                bam_handler$get_noise_ratio(chip_results$BamNames,
                                            input_results$BamNames)
            results[design_name] <- chip_results$Coverage - (input_results$Coverage * noise_ratio)
            
            # When input signal is stronger than the chip's, we'll get
            # negative values. Set the value floor to 0.
            i <- results[[design_name]] < 0
            results[[design_name]][i] <- 0
        } else {
            # If we had no input bams, return coverage as-is.
            results[design_name] <- chip_results$Coverage
        }
    }
    results
}

merge_chip = function(coverages, design) {
    result <- list()
    for (design_name in colnames(design)[-1]) {
        result[[design_name]] <- merge_reduce(coverages, design, design_name, 1)$Coverage
    }
    result
}

merge_reduce = function(coverages, design, design_name, design_value) {
    indices = design[[design_name]] == design_value
    bam_names <- design[indices, 1]
    
    list(Coverage=Reduce("+", coverages[bam_names]),
         BamNames=bam_names)
}
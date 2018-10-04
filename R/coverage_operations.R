###############################################################################
# Binning of contiguous regions (GRanges objects)                             #
###############################################################################

# Split the ranges within gr into bcount bins and returns 
# the mean coverage for each bin in a matrix where each row
# is a binned region, and each column is a bin.
get_view_means = function(gr, bcount, cov) {
    chr <- unique(as.character(GenomeInfoDb::seqnames(gr)))
    gr <- unlist(tile(gr, n=bcount))
    stopifnot(length(chr) == 1)
    views <- Views(cov[[chr]], start(gr), end(gr))
    viewMeans(views)
}


# Split coverage (an RleList) into bin_count bins over the specified 
# regions (A GRanges object). Coverages for regions on the '-' strands
# are then reversed.
bin_contiguous_regions <- function(coverage, regions, bin_count) {
  # Split regions by chromosomes, since coverages are split this way.
  grl <- split(regions, GenomeInfoDb::seqnames(regions))
  
  # Discard chromosomes where no regions are found.
  i <- vapply(grl, length, numeric(1)) > 0
  
  # Get region means on a per-chromosome basis.
  vector_means = do.call("c", lapply(grl[i], get_view_means,
                                     bcount = bin_count, cov = coverage))
  m <-  matrix(vector_means, ncol=bin_count, byrow=TRUE)

  # Reorder matrix to preserve region order.
  m <- m[GenomicRanges::findOverlaps(regions, unlist(grl), select="first"),, drop=FALSE]
  
  # Reverse rows on "-" strand
  mr <- m[,bin_count:1, drop=FALSE]
  i <-as.logical(strand(regions)=="-")
  m[i,] <- mr[i,]
  
  m
}


###############################################################################
# Binning of discontiguous regions (GRangesList objects)                      #
###############################################################################

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

# Bin one element from an RleList (x) into bin_count bins.
bin_rle_list = function(x, bin_count) { 
    viewMeans(Views(x, breakInChunks(length(x), nchunk=bin_count)))
}

# Split coverage (an RleList) into bin_count bins over the specified 
# regions (A GRangesList object).
bin_discontiguous_regions <- function(coverage, regions, bin_count) {
    rle_list = discontiguous_coverage(coverage, regions)
    bin_list = lapply(rle_list, bin_rle_list, bin_count=bin_count)
    matrix(unlist(bin_list), ncol=bin_count, byrow=TRUE)
}


###############################################################################
# Top-level binning functions to manage strand-specific coverage and          #
# branching on the type of region (contiguous GRanges, or discontiguous       #
# GRangesList                                                                 #
###############################################################################

# Call the right binning algorithm depending on the type of regions:
# contiguous (GRanges) vs discontiguous (GRangesList).
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

# Given a list of stranded coverages (+, -, *), bin them into bin_count bins over
# each element of regions.
bin_coverages_s = function(coverages, regions, bin_count) {
    # In single_strand mode, we'll only have coverage info for the undefined strand.
    if(is.null(coverages[["+"]]) && is.null(coverages[["-"]])) {
        return(bin_region_coverages(coverages[["*"]], regions, bin_count))
    } else {
        strand_split=split(regions, strand(regions))

        # Combine coverages? Star coverage should include every single read.
        # Maybe not compute star coverage and add + and - here?
        ## Figure out star coverage by combining regions.
        # star_coverage_plus = bin_region_coverages(coverages[["+"]], strand_split[["*"]], bin_count)
        # star_coverage_minus = bin_region_coverages(coverages[["-"]], strand_split[["*"]], bin_count)
        # star_coverage_star = bin_region_coverages(coverages[["*"]], strand_split[["*"]], bin_count)
        # 
        # # Recombine star coverages.
        # star_coverage = purrr::pmap(list(star_coverage_plus, star_coverage_minus, star_coverage_star), function(x, y, z) {x+y+z})
        
        coverage_list = list("+"=bin_region_coverages(coverages[["+"]], strand_split[["+"]], bin_count),
                             "-"=bin_region_coverages(coverages[["-"]], strand_split[["-"]], bin_count),
                             "*"=bin_region_coverages(coverages[["*"]], strand_split[["*"]], bin_count))
        
        recombine_regions = function(x,y,z, strand_info, bin_count) {
            results = matrix(0, nrow=length(strand_info), ncol=bin_count)
            results[strand_info=="+",] = x
            results[strand_info=="-",] = y
            results[strand_info=="*",] = z
            
            results
        }
        
        # Recombine all coverages
        purrr::pmap(coverage_list, recombine_regions, strand_info=strand(regions), bin_count=bin_count)
    }
}

###############################################################################
# Grouping of multiple coverages (bam files) into single coverages.           #
###############################################################################

# Group coverages according to design.
group_coverages_s = function(coverage_s, design, merge_operation='+', bam_handler=NULL) {
    if(is.null(coverage_s)) {
        results = NULL
    } else {
        if (merge_operation=="NCIS") {
            results <- remove_controls(coverage_s, design, bam_handler)
        } else {
            results <- merge_chip(coverage_s, design, merge_operation)
        }
    }
    
    return(results)
}

# Perform noise-reduction using the NCIS method.
remove_controls = function(coverages, design, bam_handler) {
    results <- list()
    for (design_name in colnames(design)[-1]) {
        # Add up coverage for all ChIP and all input bams.
        chip_results <- merge_reduce(coverages, design, design_name, 1, '+')
        input_results <- merge_reduce(coverages, design, design_name, 2, '+')
        
        
        if (length(input_results$BamNames) > 0) {
            # If we had input bams, perform noise reduction.
            noise_ratio <-
                bam_handler$get_noise_ratio(as.character(chip_results$BamNames),
                                            as.character(input_results$BamNames))
            results[[design_name]] <- chip_results$Coverage - (input_results$Coverage * noise_ratio)
            
            # When input signal is stronger than the chip's, we'll get
            # negative values. Set the value floor to 0.
            i <- results[[design_name]] < 0
            results[[design_name]][i] <- 0
        } else {
            # If we had no input bams, return coverage as-is.
            results[[design_name]] <- chip_results$Coverage
        }
    }
    results
}

# Merge multiple coverages within a coverage list according to the columns
# of design.
merge_chip = function(coverages, design, merge_operation) {
    result <- list()
    for (design_name in colnames(design)[-1]) {
        result[[design_name]] <- merge_reduce(coverages, design, design_name, 1, merge_operation)$Coverage
    }
    result
}

# Perform the "Reduce" operation on the elements of coverage which have
# the value design_value in column design_name of the design data-frame.
merge_reduce = function(coverages, design, design_name, design_value, merge_operation) {
    indices = design[[design_name]] == design_value
    bam_names <- design[indices, 1]
    
    reduced_coverages = Reduce("+", coverages[bam_names])
    if(merge_operation=="mean" && sum(indices) > 1) {
        reduced_coverages = reduced_coverages / sum(indices)
    }

    list(Coverage=reduced_coverages,
         BamNames=bam_names)
}

###############################################################################
# Splitting of coverages/regions based on region metadata                     #
###############################################################################

# Given a metadata data-frame and a vector of column names (split_by),
# builds a partition of rows based on all possible combinations
# of the columns identified by split_by.
split_by_metadata = function(metadata, split_by) {
    # Determine all possible values for the split_by columns.
    split_by_list = as.list(split_by)
    names(split_by_list) = split_by
    possible_values = lapply(split_by_list, function(x) {unique(metadata[[x]])})
    
    # Determine all possible combinations of the split_by column values.
    combinations = expand.grid(possible_values)

    # Split rows by iterating over value combinations.
    out_subsets=list()
    new_metadata_list = list()
    partition = rep(NA, nrow(metadata))
    for(i in 1:nrow(combinations)) {
        # Select the columns where all values match the current combination.
        selected_subset = TRUE
        for(j in 1:ncol(combinations)) {
            col_name = colnames(combinations)[j]
            col_value = combinations[i,j]
            if(!is.na(col_value)) {
                selected_subset = selected_subset & (metadata[[col_name]] == col_value) & !is.na(metadata[[col_name]])
            } else {
                selected_subset = selected_subset & is.na(metadata[[col_name]])
            }
        }
        
        # If at least one row is selected, generate a name and assign it.
        if(sum(selected_subset, na.rm=TRUE) > 0) {
            region_name = paste(combinations[i,], collapse=";")
            out_subsets[[region_name]] = selected_subset
            partition[selected_subset] = i
            
            # We'll store the combination values for later use.
            new_metadata_list[[region_name]] = combinations[i,, drop=FALSE]
        }
    }
    
    # Concatenate the new metadata.
    new_metadata=data.table::rbindlist(new_metadata_list, use.names=TRUE)
    new_metadata$split_regions = names(new_metadata_list)

    # Return the results.
    return(list(Indices=out_subsets, Metadata=new_metadata, Partition=partition))        
}

# Given a single matrix, split it according to the split_by columns of the metadata.
split_matrix = function(input_matrix, split_indices) {
    lapply(split_indices, function(indices) { input_matrix[indices,, drop=FALSE] })
}

# Given a list of matrices, split each one according to the split_by columns 
# of the metadata.
split_matrices = function(matrices, metadata, split_by) {
    split_indices = split_by_metadata(metadata, split_by)
    
    res_matrices = lapply(matrices, split_matrix, split_indices=split_indices$Indices)
    return(list(Matrices=res_matrices, Metadata=split_indices$Metadata))
}

# Given a GRange or GRangesList object, split it according to the split_by columns 
# of the metadata. GRanges are converted to GRangesList, while GRangesList
# are converted to a list of GRangesList objects.
split_regions = function(regions, metadata, split_by) {
    split_indices = split_by_metadata(metadata, split_by)

    out_regions=list()
    for(split_name in names(split_indices$Indices)) {
        out_regions[[split_name]] = regions[split_indices$Indices[[split_name]]]
    }

    # If we were splitting a GRanges, we can make this a GRangesList.
    if(is(regions, "GRanges")) {
        out_regions = GRangesList(out_regions)
    }
    
    attr(out_regions, "split_metadata") = split_indices$Metadata
    return(out_regions)
}
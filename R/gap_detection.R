# Given a set of input matrices of binned coverages,
# determines where "gaps" in coverage are found.
#
# A bin is considered to be part of a gap if:
#  - There are less than min_below regions where the coverage for
#    this bin is above threshold, and
#  - There are at least min_width adjacent bins who also fulfill
#    the first criteria.
gap_detection = function(input_matrices, threshold=0, min_above=1, min_width=1) {
    # Infer bin_count from first region/design combination.
    bin_count = ncol(input_matrices[[1]][[1]]$input)
    
    # Determine the number of single regions above the threshold for 
    # each region group/design combination.
    n_region_above=matrix(0, ncol=bin_count, nrow=0)
    for(region in names(input_matrices)) {
        for(design in names(input_matrices[[region]])) {
            n_region_above = rbind(n_region_above, apply(input_matrices[[region]][[design]]$input > threshold, 2, sum))
        }
    }
    
    # Calculate which bins have enough regions above the threshold.
    n_region_above_threshold = colSums(n_region_above)
    bin_above = n_region_above_threshold > min_above
    
    # Determine which gaps are wide enough. 
    logical_rle = Rle(bin_above)
    gap_runs = (runLength(logical_rle) >= min_width) & !runValue(logical_rle)
    
    # Reset the logical markers to TRUE where gaps weren't wide enough.
    runValue(logical_rle)[!gap_runs] = TRUE
    
    # Return gapped bins.
    return(as.logical(!logical_rle))
}
# Given a set of input matrices of binned coverages, determine
# the location of gaps and excise those columns. See gap_detection
# for details on teh gap detection algorithm.
gap_removal = function(input_matrices, threshold=0, min_above=1, min_width=1) {
    gaps = gap_detection(input_matrices, threshold, min_above, min_width)
    
    if(all(gaps)) {
        warning("Nothing left after gap detection! Gap removal will not be performed.")
        return(input_matrices)
    }
    
    # Remove the gaps
    results = list()
    for(region in names(input_matrices)) {
        results[[region]] = list()
        for(design in names(input_matrices[[region]])) {
            results[[region]][[design]] = list(input=input_matrices[[region]][[design]]$input[,!gaps, drop=FALSE])
        }
    }

    return(results)
}
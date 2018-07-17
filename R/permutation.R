# Given a vector x, resamples it sample_count time and returns
# the mean and confidence intervals at the alpha level.
calculate_bin_ci = function(x, sample_count, alpha, sample_indices=NULL) { 
    if(is.null(sample_indices)) {
        sample_indices = sample.int(length(x), size=length(x)*sample_count, replace=TRUE)
    }
    
    # Resample into a matrix of length(x) rows and sample_count columns.
    sampled = matrix(x[sample_indices], ncol=sample_count);
    
    # Calculate the column means, which are the means of each resampling.
    means = colMeans(sampled);
    
    # Put the results in a named vector.
    return(c(mean(x), quantile(means, c(alpha/2, 1-(alpha/2)))))
}

# Given a list with elements Region, Design and Matrix, resamples all columns
# of Matrix sample_count times and calculate confidence intervals of the means at level alpha.
# The results are stored as a data-frame with the additional design and region columns.
calculate_matrix_ci = function(x, sample_count, alpha, reuse) {
    if(reuse) {
        sample_indices = sample.int(nrow(x$Matrix), size=nrow(x$Matrix)*sample_count, replace=TRUE)
    } else {
        sample_indices = NULL
    }
    
    # Resample and calculate CIs for all columns of the matrix.
    res = t(apply(x$Matrix, 2, calculate_bin_ci, sample_count=sample_count, alpha=alpha, sample_indices=sample_indices))
    
    # Format the resulting data-frame correctly.
    colnames(res) = c("value", "qinf", "qsup")
    res = data.frame(res)
    res$design = x$Design
    res$region = x$Region
    
    res
}

calculate_matrices_ci = function(matrices, sample_count, alpha, resampling_mode, parallel_job=NULL) {
    # Get coverage matrices, and reformat them into a flat list
    # so each matrix can be processed in parralel.
    matrix_list = list()
    i=1
    for(design in names(matrices)) {
        for(region in names(matrices[[design]])) {
            matrix_list[[i]] = list(Region=region, Design=design, Matrix=matrices[[design]][[region]])
            i = i + 1
        }
    }
    
    # Perform resampling in parralel.
    if(resampling_mode=="profile") {
        reuse=TRUE
    } else if(resampling_mode=="bin") {
        reuse=FALSE
    } else {
        stop("resampling_mode must be 'profile' or 'bin'")
    }
    
    if(!is.null(parallel_job)) {
        ci <- parallel_job$launch_job(
                    data = matrix_list,
                    FUN = calculate_matrix_ci,
                    sample_count = sample_count,
                    alpha = alpha,
                    reuse=reuse)
    } else {
        ci = lapply(matrix_list, calculate_matrix_ci, sample_count = sample_count, alpha = alpha, reuse=reuse)
    }
    
    # Concatenate resampling results and add bin column.
    res = data.table::rbindlist(ci, idcol=NULL, use.names=TRUE, fill=FALSE)
    res$bin = 1:ncol(matrix_list[[1]]$Matrix)
    
    res
}

add_metadata_to_ci = function(ci_df, region_metadata, design_metadata) {
     # Add metadata.
     results <- as.data.frame(ci_df)
     results$group <- as.factor(paste(results$design, results$region, sep="_"))
     results = dplyr::left_join(results, region_metadata, by=c(region="split_regions"))      
     results = dplyr::left_join(results, design_metadata, by=c(design="design"))      
     
     return(results)   
}
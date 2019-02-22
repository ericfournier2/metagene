validate_design_format = function(design) {
    if(!is.data.frame(design)) {
        stop("design must be a data.frame object, NULL or NA")
    }

    # Validate that we have enough columns and that they are of the right types.
    if(ncol(design) < 2) {
        stop("design must have at least 2 columns")
    }
    if (!(is.character(design[,1]) || is.factor(design[,1]))) {
        stop("The first column of design must be BAM filenames")
    }
    if (!all(apply(design[, -1, drop=FALSE], MARGIN=2, is.numeric))) {
        stop(paste0("All design column, except the first one,",
                        " must be in numeric format"))
    }
}

validate_alpha = function(alpha) {
    stopifnot(is.numeric(alpha))
    stopifnot(alpha >= 0 & alpha <= 1)
}

validate_bin_count = function(bin_count) {
    if (!is.null(bin_count)) {
        if (!is.numeric(bin_count) || bin_count <= 0 || as.integer(bin_count) != bin_count) {
            stop("bin_count must be a positive integer")
        }
    }
}

validate_sample_count = function(sample_count) {
    stopifnot(is.numeric(sample_count))
    stopifnot(sample_count >= 0)
    stopifnot(as.integer(sample_count) == sample_count)
}

validate_avoid_gaps = function(avoid_gaps) {
    stopifnot(is.logical(avoid_gaps))
}

validate_gaps_threshold = function(gaps_threshold) {
    stopifnot(is.numeric(gaps_threshold))
    stopifnot(gaps_threshold >= 0)
}

validate_normalization = function(normalization) {
    if (!is.null(normalization) && normalization != "RPM" && normalization != "NCIS") {
        stop('normalization must be NULL, "RPM" or "NCIS".')
    }
}

validate_flip_regions = function(flip_regions) {
    if (!is.logical(flip_regions)) {
        stop("flip_regions must be a logical.")
    }
}

validate_assay = function(assay) {
    # Check parameters validity
    if (!is.character(assay)) {
        stop("verbose must be a character value")
    }
    assayTypeAuthorized <- c('chipseq', 'rnaseq')
    if (!(tolower(assay) %in% assayTypeAuthorized)) {
        stop("assay values must be one of 'chipseq' or 'rnaseq'")
    }
}

validate_verbose = function(verbose) {
    if (!is.logical(verbose)) {
        stop("verbose must be a logicial value (TRUE or FALSE)")
    }       
}

validate_force_seqlevels = function(force_seqlevels) {
    if (!is.logical(force_seqlevels)) {
        stop(paste("force_seqlevels must be a logicial ",
                    "value (TRUE or FALSE)",sep=""))
    }        
}

validate_padding_size = function(padding_size) {
    if (!(is.numeric(padding_size) || is.integer(padding_size)) ||
        padding_size < 0 || as.integer(padding_size) != padding_size) {
        stop("padding_size must be a non-negative integer")
    }       
}

validate_cores = function(cores) {
    isBiocParallel = is(cores, "BiocParallelParam")
    isInteger = ((is.numeric(cores) || is.integer(cores)) &&
                    cores > 0 && as.integer(cores) == cores)
    if (!isBiocParallel && !isInteger) {
        stop("cores must be a positive integer or a BiocParallelParam instance.")
    }        
}

validate_bam_files = function(bam_files) {
    validate_bam_files_format(bam_files)
    validate_bam_files_values(bam_files)
}

validate_bam_files_format = function(bam_files) {
    if (!is.vector(bam_files, "character")) {
        stop("bam_files must be a vector of BAM filenames.")
    }      
}

validate_bam_files_values = function(bam_files) {
    if (!all(sapply(bam_files, file.exists))) {
        stop("At least one BAM file does not exist.")
    }        
}   
     
validate_combination = function(params) {
    if(params$assay=="chipseq" && is.null(params$bin_count)) {
        stop("bin_count cannot be NULL in chipseq assays.")
    }
}

validate_region_mode = function(region_mode) {
    if(!(region_mode %in% c("auto", "separate", "stitch"))) {
        stop("region_mode must be 'auto', 'separate' or 'stitch'")
    }
}

validate_extend_reads = function(extend_reads) {
    stopifnot(is.numeric(extend_reads))
    stopifnot(extend_reads >= 0)
}
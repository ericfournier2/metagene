# Class used to manage parallel jobs
Parallel_Job <- R6Class("Parallel_Job",
    public = list(
        parameters = list(),
        BPPARAM = SerialParam(),
        initialize = function(cores = SerialParam()) {
            self$set_core_count(cores)
        },
        launch_job = function(data, FUN, ...) {
            res <- BiocParallel:::bplapply(data, FUN, BPPARAM = self$BPPARAM,
                                           ...)
            # Check for errors
            if (!all(bpok(res))) {
                i <- grepl("stop", attr(res[[1]], "traceback"))
                i <- which(i)[1] # We return first error message
                msg <- strsplit(attr(res[[1]], "traceback")[i], "[()]")[[1]][2]
                msg <- substring(msg, 2, nchar(msg)-1)
                stop(msg)
            }
            res
        },
        get_core_count = function() {
            cores <- self$parameters[["cores"]]
            if(is.numeric(cores)) {
                cores
            } else {
                bpworkers(cores)
            }
        },
        set_core_count = function(cores) {
            # Stop any previous thread pool. I can't really fathom
            # a situation where people just keep setting new BPPARAMs
            # until we run out of threads or memory, but better safe than sorry.
            tryCatch({
                bpstop(self$BPPARAM)
            })            
            
            # Note: cores can be numeric or BiocParallelParam instance
            #       BPPARAM is always a BiocParallelParam instance
            self$BPPARAM <- private$set_bpparam(cores)
            self$parameters[["cores"]] <- cores            
        }
    ),
    private = list(
        set_bpparam = function(cores) {
			if (is.numeric(cores)) {
                # The number of cores has to be a positive integer
                if(as.integer(cores) != cores || cores <= 0) {
                    msg <- "cores must be positive numeric or "
                    msg <- paste0(msg, "BiocParallelParam instance.")
                    stop(msg)
                }
                if (cores == 1) {
                    BPPARAM <- SerialParam()
                } else if (.Platform$OS.type != "unix") {
                    BPPARAM <- SnowParam(workers = cores)
                } else {
                    BPPARAM <- MulticoreParam(workers = cores)
                }
            } else {
            # Must be one of the BiocParallelParam class
                if (!is(cores, "BiocParallelParam")) {
                    msg <- "cores must be positive numeric or "
                    msg <- paste0(msg, "BiocParallelParam instance.")
                    stop(msg)
                }
                BPPARAM <- cores
            }
            BPPARAM
        }
    )
)

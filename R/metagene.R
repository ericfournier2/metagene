#' A class to manage metagene analysis.
#'
#' This class will allow to load, convert and normalize alignments and regions
#' files/data. Once the data is ready, the user can then choose to produce
#' metagene plots on the data or some subset of it.
#'
#' @section Constructor:
#' \describe{
#'    \item{}{\code{mg <- metagene$new(regions, bam_files, padding_size = 0,
#'                            cores = SerialParam(), verbose = FALSE,
#'                            force_seqlevels = FALSE, paired_end = FALSE,
#'                            assay = 'chipseq'))}}
#'    \item{regions}{Either a \code{vector} of filenames, a \code{GRanges} object 
#'                    or a \code{GRangesList} object. Supported file formats are
#'                    BED, narrowPeak, broadPeak, gff and gtf. For rnaseq assays,
#'                    a GRanges object represents a single gene, while multiple
#'                    genes are represented by a \code{GRangesList} object, where
#'                    each individual \code{GRanges} is a set of exons.}
#'    \item{bam_files}{A \code{vector} of BAM filenames. The BAM files must be
#'                    indexed. i.e.: if a file is named file.bam, there must
#'                    be a file named file.bam.bai or file.bai in the same 
#'                    directory.}
#'    \item{padding_size}{The regions will be extended on each side by the
#'                        value of this parameter. The padding_size must be a
#'                        non-negative integer. Default = 0.}
#'    \item{cores}{The number of cores available to parallelize the analysis.
#'                Either a positive integer or a \code{BiocParallelParam}.
#'                Default: \code{SerialParam()}.}
#'    \item{verbose}{Print progression of the analysis. A logical constant.
#'                    Default: \code{FALSE}.}
#'    \item{force_seqlevels}{If \code{TRUE}, Remove regions that are not found
#'                in bam file header. Default: \code{FALSE}. TRUE and FALSE
#'                respectively correspond to pruning.mode = "coarse"
#'                and "error" in ?seqinfo.}
#'    \item{paired_end}{If \code{TRUE}, metagene will deal with paired-ended 
#'                data. If \code{FALSE}, single-ended data are expected. 
#'                Default: \code{FALSE}}
#'    \item{assay}{\code{'chipseq'} or \code{'rnaseq'} Default: \code{'chipseq'}}
#'    \item{strand_specific}{If \code{TRUE}, only reads which align to the same 
#'                           strand as those specified in \code{regions} will
#'                           count toward coverage for that region. Useful for RNA-seq
#'                           profiles generated from strand-specific libraries, such
#'                           as Illumina TruSeq. Default: \code{'FALSE'}}
#'    \item{paired_end_strand_mode}{\code{'1'} or \code{'2'}. In paired-end mode
#'                                  indicates which read in a pair sets the pair's strand.
#'                                  If \code{1}, this is the first read (This should be used
#'                                  with directional protocols such as Directional Illumina 
#'                                  (Ligation) or Standard SOLiD).
#'                                  If \code{2}, this is the second read (This should be used
#'                                  with directional protocols such as dUTP, NSR, NNSR, 
#'                                  or Illumina stranded TruSeq PE).
#'                                  Ignored if either paired_end or strand_specific is FALSE.
#'                                  Default: \code{'2'}}
#' }
#'
#'    \code{metagene$new} returns a \code{metagene} object that contains the
#'        coverages for every BAM files in the regions from the \code{regions}
#'        param.
#'
#' @return
#' \code{metagene$new} returns a \code{metagene} object which contains the
#' normalized coverage values for every regions in all specified BAM files.
#'
#' @section Methods:
#' \describe{
#'    \item{}{\code{mg$plot(region_names = NULL, design_names = NULL,
#'                title = NULL, x_label = NULL)}}
#'    \item{region_names}{The names of the regions to extract. If \code{NULL},
#'                        all the regions are returned. Default: \code{NULL}.}
#'    \item{design_names}{The names of the experiments to extract. If a design
#'            was added to the \code{metagene} object, \code{design_names}
#'            correspond to the column names in the design, otherwise
#'            \code{design_names} corresponds to the BAM name or the BAM
#'            filename. If \code{NULL}, all the experiments are
#'            returned. Default: \code{NULL}.}
#'    \item{title}{A title to add to the graph. If \code{NULL}, will be
#'                        automatically created. Default: NULL}
#'    \item{x_label}{X-axis label to add to the metagene plot. If \code{NULL},
#'                    metagene will use generic label. Default: \code{NULL}.}
#' }
#' \describe{
#'    \item{}{\code{mg$produce_table(design, bin_count, noise_removal,
#'                normalization, flip_regions, bin_size = NULL}}
#'    \item{design}{A \code{data.frame} that describe the experiment to plot.
#'            see \code{plot} function for more details. \code{NA} can 
#'            be used to keep the previous design value. Default: \code{NA}.}
#'    \item{bin_count}{The number of bins to create. \code{NA} can be used to
#'                    keep the previous bin_count value. For ChIP-Seq analyses,
#'                    a bin_count of 100 will be used if no value is specified.
#'                    For RNA-seq analyses, no binning will occur unless a bin_count
#'                    is explicitly specified. Default: \code{NA}.}
#'    \item{noise_removal}{The algorithm to use to remove control(s). Possible
#'                        values are \code{NA}, \code{NULL} or "NCIS". By
#'                        default, value is \code{NULL}. Use \code{NA} to keep
#'                        the previous \code{noise_removal} value (i.e. if
#'                        \code{produce_table} was called before). See
#'                        Liand and Keles 2012 for the NCIS algorithm.}
#'    \item{normalization}{The algorithm to use to normalize samples. Possible
#'                        values are \code{NA}, \code{NULL} or "RPM". By
#'                        default, value is \code{NULL} and no normalization
#'                        will be performed. Use \code{NA} keep
#'                        previous \code{normalization} value (i.e. if
#'                        \code{produce_table} was called before).}
#'    \item{flip_regions}{Should regions on the negative strand be flipped?
#'                        Default: \code{FALSE}.}
#'    \item{bin_size}{Deprecated.}
#' }
#' \describe{
#'    \item{}{\code{mg$produce_data_frame(alpha = 0.05, sample_count = 1000, 
#'                                avoid_gaps = FALSE, gaps_threshold = 0)}}
#'    \item{alpha}{The range of the estimation to be shown with the ribbon.
#'                \code{1 - alpha / 2} and \code{alpha / 2} will be used.
#'                Default: 0.05.}
#'    \item{sample_count}{The number of draws to perform in the bootstrap
#'                        calculations. Default: 1000.}
#'    \item{avoid_gaps}{Removes regions below the coverage threshold specified 
#'                      by \code{gaps_threshold}, and refits the data_frame 
#'                      accordingly. Default: \code{FALSE}.}
#'    \item{gaps_threshold}{Sets the threshold used for the gap removal procedure.
#'                          Regions with values values <= gaps_threshold will be removed.
#'                          Ignored if avoid_gaps is \code{FALSE}. Default: 0.}
#' }
#' \describe{
#'    \item{}{mg$get_params()}
#' }
#' \describe{
#'    \item{}{mg$get_design()}
#' }
#' \describe{
#'    \item{}{mg$get_regions(region_names = NULL)}
#'    \item{region_names}{The names of the regions to extract. If \code{NULL},
#'                        all the regions are returned. Default: \code{NULL}.}
#' }
#' \describe{
#'    \item{}{mg$get_table = function()}
#' }
#' \describe{
#'    \item{}{mg$get_matrices = function()}
#' }
#' \describe{
#'    \item{}{mg$get_data_frame(region_names = NULL, design_names = NULL)}
#'    \item{region_names}{The names of the regions to extract. If \code{NULL},
#'                        all the regions are returned. Default: \code{NULL}.}
#'    \item{design_names}{The names of the experiments to extract. If a design
#'            was added to the \code{metagene} object, \code{design_names}
#'            correspond to the column names in the design, otherwise
#'            \code{design_names} corresponds to the BAM name or the BAM
#'            filename. If \code{NULL}, all the experiments are
#'            returned. Default: \code{NULL}.}
#' }
#' \describe{
#'    \item{}{get_plot = function()}
#' }
#' \describe{
#'    \item{}{get_raw_coverages = function(filenames)}
#'    \item{filenames}{The name of the files from which raw coverages should
#'                     be extracted. Can be the filenames with extensions or
#'                     the bam names, if named bam files were used during the
#'                     creation of the metagene object). If \code{NULL}, 
#'                     returns the coverages for all bam files. Default: 
#'                     \code{NULL}.}
#' }
#' \describe{
#'    \item{}{get_normalized_coverages = function(filenames)}
#'    \item{filenames}{The name of the files from which normalized coverages
#'                     (in RPM)  should be extracted. Can be the filenames with
#'                     extensions or the bam names, if named bam files were 
#'                     used during the creation of the metagene object). If 
#'                     \code{NULL}, returns the coverages for all bam files. 
#'                     Default: \code{NULL}.}
#' }
#' \describe{
#'    \item{}{\code{mg$export(bam_file, region, file)}}
#'    \item{bam_file}{The name of the bam file to export.}
#'    \item{region}{The name of the region to export.}
#'    \item{file}{The name of the ouput file.}
#' }
#' \describe{
#'    \item{}{\code{mg$add_design(design = NULL, check_bam_files = FALSE)}}
#'    \item{design}{A \code{data.frame} that describes the experiment to plot.
#'            See \code{plot} function for more details. \code{NA} can be
#'            used to keep the previous design value. Default: \code{NA}.}
#'    \item{check_bam_files}{Force check that all bam files from the first
#'                           column of the design are present in the current
#'                           metagene object. Default: \code{FALSE}}
#' }
#' \describe{
#'    \item{}{\code{mg$unflip_regions()}}
#' }
#'
#' \describe{
#'    \item{}{\code{mg$flip_regions()}}
#' }
#' @examples
#' region <- get_demo_regions()[1]
#' bam_file <- get_demo_bam_files()[1]
#' mg <- metagene$new(regions = region, bam_files = bam_file)
#' \dontrun{
#'    df <- metagene$plot()
#' }
#'
#' @importFrom R6 R6Class
#' @importFrom data.table data.table
#' @export
#' @format A metagene experiment manager

metagene <- R6Class("metagene",
    public = list(
    # Methods
        initialize = function(regions, bam_files, padding_size = 0,
                                cores = SerialParam(), verbose = FALSE,
                                force_seqlevels = FALSE, paired_end = FALSE,
                                assay = 'chipseq', strand_specific=FALSE,
                                paired_end_strand_mode=2) {

            # Validate the format of bam_files, since it is used to preprocess certain
            # parameters before initialization.
            private$validate_bam_files_format(bam_files)
            
            # Initialize parameter handler.
            private$ph <- parameter_manager$new(
                param_values=list(
                    design=private$get_complete_design(bam_files),
                    bam_files=private$name_from_path(bam_files),
                    padding_size=padding_size,
                    verbose=verbose,
                    force_seqlevels=force_seqlevels,
                    paired_end=paired_end,
                    assay=tolower(assay),
                    strand_specific=strand_specific,
                    paired_end_strand_mode=paired_end_strand_mode,
                    normalization=NULL,
                    noise_removal=NULL,
                    avoid_gaps=FALSE,
                    gaps_threshold=0,
                    bin_count=100,
                    flip_regions=FALSE,
                    alpha=0.05,
                    sample_count=1000,
                    bam_name=NULL),
                param_validations=list(
                    design=private$validate_design,
                    bam_files=private$validate_bam_files,
                    padding_size=private$validate_padding_size,
                    verbose=private$validate_verbose,
                    force_seqlevels=private$validate_force_seqlevels,
                    paired_end=private$validate_paired_end,
                    assay=private$validate_assay,
                    strand_specific=private$validate_strand_specific,
                    paired_end_strand_mode=private$validate_paired_end_strand_mode,
                    normalization=private$validate_normalization,
                    noise_removal=private$validate_noise_removal,
                    avoid_gaps=private$validate_avoid_gaps,
                    gaps_threshold=private$validate_gaps_threshold,
                    bin_count=private$validate_bin_count,
                    flip_regions=private$validate_flip_regions,
                    alpha=private$validate_alpha,
                    sample_count=private$validate_sample_count),
                overall_validation=private$validate_combination)
                
            # Prepare objects for parralel processing.
            private$validate_cores(cores)
            private$parallel_job <- Parallel_Job$new(cores)
            
            # Prepare bam files
            private$print_verbose("Prepare bam files...")
            private$bam_handler <- Bam_Handler$new(private$ph$get("bam_files") , cores = cores,
                                        paired_end = paired_end,
                                        strand_specific=strand_specific,
                                        paired_end_strand_mode=paired_end_strand_mode)

            # Prepare regions
            private$print_verbose("Prepare regions...")
            private$regions <- private$prepare_regions(regions, private$ph$get("assay"))

            # Parse bam files
            private$print_verbose("Calculating coverages...\n")
            private$coverages <- private$produce_coverages()    
        },
        get_bam_count = function(filename) {
            # Parameters validation are done by Bam_Handler object
            private$bam_handler$get_aligned_count(filename)
        },
        get_params = function() {
            private$ph$get_all()
        },
        get_design = function() {
            private$ph$get("design")
        },
        get_regions = function(region_names = NULL) {
            if (is.null(region_names)) {
                private$regions
            } else {
                new_names <- tools::file_path_sans_ext(basename(region_names))
                region_names <- new_names
                stopifnot(all(region_names %in% names(private$regions)))
                private$regions[region_names]
            }
        },
        get_table = function() {
            if (length(private$table) == 0) { 
                return(NULL)
            }
            assay = private$ph$get("assay")
            if (assay == 'chipseq'){
                return(copy(private$table))
            } else if (assay == 'rnaseq'){
                if(is.null(private$ph$get('bin_count'))) {
                    return(copy(private$table[,-"bin"]))
                } else {
                    return(copy(private$table))
                }
            }
        },
        get_matrices = function() {
            if (is.null(self$get_table())){
                return(NULL)
            }
            assay = private$ph$get("assay")
            if (assay == 'chipseq') {
                matrices <- list()
                nbcol <- private$ph$get("bin_count")
                nbrow <- vapply(self$get_regions(), length, numeric(1))
                for (regions in names(self$get_regions())) {
                    matrices[[regions]] <- list()
                    for (design_name in colnames(self$get_design())[-1]) {
                        matrices[[regions]][[design_name]] <- list()
                        matrices[[regions]][[design_name]][["input"]] <- 
                                matrix(private$table[region == regions & 
                                design == design_name,]$value, 
                                nrow=nbrow[regions], ncol=nbcol, byrow=TRUE)
                    }
                }
                return (matrices)
            } else {
                stop(paste('unsupported function for assay of type',
                        assay,
                        '. Only available for "chipseq" assay.')) 
            }
        },
        get_data_frame = function(region_names = NULL, design_names = NULL) {
            if (nrow(private$df) == 0) {
                NULL
            } else if (is.null(region_names) && is.null(design_names)) {
                return(copy(private$df))
            } else {
                if (!is.null(region_names)) {
                    stopifnot(is.character(region_names))
                    stopifnot(all(region_names %in% unique(private$table$region)))
                } else {
                    region_names <- names(private$regions)
                }
                if (!is.null(design_names)) {
                    stopifnot(is.character(design_names))
                    stopifnot(all(design_names %in% unique(private$table$design)))
                } else {
                    design_names <- colnames(private$ph$get("design"))[-1]
                }
                i <- (private$df$region %in% region_names &
                      private$df$design %in% design_names)
                return(copy(private$df[i,]))
            }
        },
        get_plot = function() {
            if (is.character(private$graph)) {
                NULL
            } else {
                private$graph
            }
        },
        get_raw_coverages = function(filenames = NULL) {
            if(!private$ph$get('strand_specific')) {
                return(private$get_raw_coverages_internal(filenames)[['*']])
            } else {
                return(private$get_raw_coverages_internal(filenames))
            }        
        },
        get_normalized_coverages = function(filenames = NULL) {
            if(!private$ph$get('strand_specific')) {
                return(private$get_normalized_coverages_internal(filenames)[['*']])
            } else {
                return(private$get_normalized_coverages_internal(filenames))
            }
        },
        #get_design_coverages = function(design=NA, noise_removal=NA, normalization=NA) {
        #    if(design_coverage_need_update(design, normalization, noise_removal)) {
        #        # Get the correct parameters.
        #        design <- private$fetch_design(design)
        #        noise_removal = private$get_param_value(noise_removal, "noise_removal")
        #        normalization <- private$get_param_value(normalization, "normalization")
        #        
        #        # Get raw coverages
        #        coverages <- private$coverages
        #        
        #        # Normalize if required.
        #        if (!is.null(normalization)) {
        #            coverages <- private$get_normalized_coverages(coverages)
        #            message('Normalization done')
        #        }
        #        
        #        # Merge the various samples, removing noise if it was requested.
        #        if (!is.null(noise_removal)) {
        #            coverages <- private$remove_controls(coverages, design)
        #        } else {
        #            coverages <- private$merge_chip(coverages, design)
        #        }
        #    
        #        private$design_coverages <- coverages
        #    }
        #    return(private$design_coverages)
        #},
        add_design = function(design, check_bam_files = FALSE) {
            design = private$clean_design(design, private$ph$get("bam_files"))
            
            # Update the design and remove invalidated table and data-frame.
            if(private$ph$update_params(design)) {
                private$table <- NULL
                private$df <- NULL
            }
        },
        produce_table = function(design = NA, bin_count = NA, bin_size = NULL,
                                noise_removal = NA, normalization = NA,
                                flip_regions = FALSE) {
            if (!is.null(bin_size)) {
                warning("bin_size is now deprecated. Please use bin_count.")
            }

            private$validate_flip_regions(flip_regions)
            design = private$clean_design(design, private$ph$get("bam_files"))
            
            if(private$ph$update_params(design, bin_count, noise_removal, normalization)) {
                private$table = NULL
            }
            
            if (is.null(private$table)) {            
                # Normalize if necessary.
                if (!is.null(private$ph$get("normalization"))) {
                    coverages <- private$get_normalized_coverages_internal()
                    message('Normalization done')
                } else {
                    coverages <- private$get_raw_coverages_internal()
                }

                # Loop over all strands, building a table for each.
                table_list = list()
                for(strand_name in c('+', '-', '*')) {
                    table_list[[strand_name]] = private$produce_strand_table(coverages[[strand_name]],
                                                                             private$ph$get("assay"), 
                                                                             private$ph$get("design"), 
                                                                             self$get_regions(), 
                                                                             private$ph$get("noise_removal"), 
                                                                             private$ph$get("bin_count"))
                }

                # Merge coverage tables.
                private$table = do.call(rbind, table_list)
                private$df <- NULL
            }

            if (flip_regions) {
                self$flip_regions()
            } else {
                self$unflip_regions()
            }
            invisible(self)
        },
        produce_data_frame = function(alpha = 0.05, sample_count = 1000, 
                                                    avoid_gaps = FALSE, 
                                                    bam_name = NULL, 
                                                    gaps_threshold = 0) {
            #arguments checking
            if (!is.null(bam_name)){
                stopifnot(is.character(bam_name))
                bam_name <- tools::file_path_sans_ext(basename(bam_name))
                bam_names <- tools::file_path_sans_ext(basename(
                                                private$ph$get("bam_files")))
                if (!bam_name %in% bam_names){
                    stop(paste("bam_name argument is not one of bam_names",
                                        "provided to the metagene object"))
                }
            }

            # If parameters hae been updated, set the table to NULL.
            if (private$ph$update_params(alpha, sample_count, avoid_gaps, bam_name, gaps_threshold)) {
                private$df <- NULL
            }
            
            if (is.null(private$df)) {
                # 1. Get the correctly formatted table
                if (is.null(self$get_table())) {
                    self$produce_table()
                }
            
                # 2. Produce the data.frame 
                private$df <- private$produce_data_frame_internal(input_table=self$get_table(),
                    alpha=alpha, sample_count=sample_count, avoid_gaps=avoid_gaps, 
                    gaps_threshold=gaps_threshold, bam_name=bam_name,
                    assay=private$ph$get('assay'), input_design=self$get_design,  
                    bin_count=private$ph$get('bin_count'))
                    
                invisible(self)
            }
        },
        plot = function(region_names = NULL, design_names = NULL, title = NULL,
                        x_label = NULL) {
            # 1. Get the correctly formatted table
            if (length(private$table) == 0) {
                self$produce_table()
            }

            # 2. Produce the data frame
            if (nrow(private$df) == 0) {
                self$produce_data_frame()
            }
            df <- self$get_data_frame(region_names = region_names,
                                    design_names = design_names)
            # 3. Produce the graph
            if (is.null(title)) {
                title <- paste(unique(private$df[["group"]]), collapse=" vs ")
            }
            p <- private$plot_graphic(df = df, title = title, 
                                        x_label = x_label)
            print(p)
            private$graph <- p
            invisible(self)
        },
        export = function(bam_file, region, file) {
            region <- tools::file_path_sans_ext(basename(region))
            region <- private$regions[[region]]
            param <- Rsamtools::ScanBamParam(which = region)
            alignments <- GenomicAlignments::readGAlignments(bam_file,
                                                                param = param)
            weight <- 1 - private$bam_handler$get_rpm_coefficient(bam_file)
            seqlevels(alignments) <- seqlevels(region)
            # TODO: don't use the weight param of coverage
            coverage <- GenomicAlignments::coverage(alignments, weight=weight)
            rtracklayer::export(coverage, file, "BED")
            invisible(coverage)
        },
        flip_regions = function() {
            if (!private$ph$get("flip_regions")) {
                private$flip_table()
                private$ph$set("flip_regions", TRUE)
            }
            invisible(self)
        },
        unflip_regions = function() {
            if (private$ph$get("flip_regions")) {
                private$flip_table()
                private$ph$set("flip_regions", FALSE)
            }
            invisible(self)
        }
    ),
    private = list(
        params = list(),
        regions = GRangesList(),
        table = NULL,
        design = data.frame(),
        coverages = list(),
        design_coverages = list(),
        df = data.frame(),
        graph = "",
        bam_handler = "",
        parallel_job = "",
        ph=NULL,
        
        #design_coverage_need_update = function(design, normalization, noise_removal) {
        #    return((length(private$design_coverages)==0) ||
        #           private$have_params_changed(design, normalization, noise_removal))
        #},
        print_verbose = function(to_print) {
            if (private$ph$get("verbose")) {
                cat(paste0(to_print, "\n"))
            }
        },
        get_subtable = function(coverages, region, bcount) {
            gr <- private$regions[[region]]
            grl <- split(gr, GenomeInfoDb::seqnames(gr))
            i <- vapply(grl, length, numeric(1)) > 0
            do.call("c", lapply(grl[i], private$get_view_means,
                                bcount = bcount, cov = coverages))
        },
        get_view_means = function(gr, bcount, cov) {
            chr <- unique(as.character(GenomeInfoDb::seqnames(gr)))
            gr <- intoNbins(gr, bcount)
            stopifnot(length(chr) == 1)
            views <- Views(cov[[chr]], start(gr), end(gr))
            viewMeans(views)
        },
        prepare_regions = function(regions, assay) {
            if (class(regions) == "character") {
                # Validation specific to regions as a vector
                if (!all(sapply(regions, file.exists))) {
                    stop("regions is a list of files, but some of those files do not exist.")
                }
                regions = private$import_regions_from_disk(regions)
            } else if (class(regions) == "GRanges") {
                regions <- GRangesList("regions" = regions)
            } else if (class(regions) == "list") {
                regions <- GRangesList(regions)
            } else if (!is(regions, "GRangesList")) {
                stop(paste0("regions must be either a vector of BED ",
                            "filenames, a GRanges object or a GrangesList object"))            
            }
            
            if (is.null(names(regions))) {
                names(regions) <- sapply(seq_along(regions), function(x) {
                    paste("region", x, sep = "_")
                    })        
            }
            
            if(assay=="rnaseq" && is(regions, "GRanges")) {
                if(length(unique(seqnames(regions))) > 1) {
                    stop(paste0("For rnaseq assays, regions should be a ",
                                "GRangesList of transcripts, or a GRanges ",
                                " object representing a single transcript. ",
                                "Here regions spans several seqnames, indicating ",
                                "it might include many transcripts."))
                }
            }

            
            if (private$ph$get('assay') == "rnaseq"){
                stopifnot(all(sum(width(GenomicRanges::reduce(private$regions)))
                            == sum(width(private$regions))))
            }
            # TODO: Check if there is a id column in the mcols of every ranges.
            #    If not, add one by merging seqnames, start and end.

            GRangesList(lapply(regions, function(x) {
                # Add padding
                padding_size = private$ph$get("padding_size")
                start(x) <- start(x) - padding_size
                start(x)[start(x) < 0] <- 1
                end(x) <- end(x) + padding_size
                # Clean seqlevels
                x <- sortSeqlevels(x)
                #seqlevels(x) <- unique(as.character(seqnames(x)))
                x
            }))
        },
        produce_coverages = function() {
            regions <- GenomicRanges::reduce(BiocGenerics::unlist(private$regions))
            bam_files = private$ph$get("bam_files")
            res <- private$parallel_job$launch_job(
                        data = bam_files,
                        FUN = private$bam_handler$get_coverage,
                        regions = regions,
                        force_seqlevels= private$ph$get("force_seqlevels"),
                        simplify=FALSE)

            names(res) <- names(bam_files)
            
            # Turn res inside out so that strand is at the top level,
            # and bam files on the second.
            res = list('+'=purrr::map(res, '+'),
                       '-'=purrr::map(res, '-'),
                       '*'=purrr::map(res, '*'))
            replace_nulls = function(x) { 
                if(all(purrr::map_lgl(x, is.null))) {
                    return(NULL)
                } else {
                    return(x)
                }
            }
            res = lapply(res, replace_nulls)

                       
            sortseq_or_null <- function(x) {
                if(is.null(x)) {
                    return(x)
                } else {
                    return(lapply(x, GenomeInfoDb::sortSeqlevels))
                }
            }
            
            lapply(res, sortseq_or_null)
        },
        plot_graphic = function(df, title, x_label) {
            # Prepare x label
            assay = private$ph$get("assay")
            if (is.null(x_label)) {
                if (assay == "chipseq") {
                    x_label <- "Distance in bins"
                } else if (assay == "rnaseq") {
                    x_label = ifelse(is.null(private$ph$get("bin_count")),
                                     "Distance in nucleotides",
                                     "Distance in bins")
                }
            }

            # Prepare y label
            y_label <- "Mean coverage"
            if (is.null(private$ph$get("normalization"))) {
                y_label <- paste(y_label, "(raw)")
            } else {
                y_label <- paste(y_label, "(RPM)")
            }
            
            # Produce plot
            p <- plot_metagene(df) +
                ylab(y_label) +
                xlab(x_label) +
                ggtitle(title)
            p
        },
        remove_controls = function(coverages, design) {
            results <- list()
            for (design_name in colnames(design)[-1]) {
                # Add up coverage for all ChIP and all input bams.
                chip_results <- private$merge_reduce(coverages, design, design_name, 1)
                input_results <- private$merge_reduce(coverages, design, design_name, 2)
                
                
                if (length(input_results$BamNames) > 0) {
                    # If we had input bams, perform noise reduction.
                    noise_ratio <-
                        private$bam_handler$get_noise_ratio(chip_results$BamNames,
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
        },
        merge_chip = function(coverages, design) {
            result <- list()
            for (design_name in colnames(design)[-1]) {
                result[[design_name]] <- private$merge_reduce(coverages, design, design_name, 1)$Coverage
            }
            result
        },
        merge_reduce = function(coverages, design, design_name, design_value) {
            indices = design[[design_name]] == design_value
            bam_files <- as.character(design[,1][indices])
            bam_names <- private$get_bam_names(bam_files)
            
            list(Coverage=Reduce("+", coverages[bam_names]),
                 BamNames=bam_names)
        },
        flip_table = function() {
            assay = private$ph$get('assay')
            if(!all(private$table[,length(levels(as.factor(strand))), 
                                    by=region][,2] == 1) &
                                    assay == 'rnaseq'){
                stop(paste('Strands of exons in one gene/region',
                                'must have the same sign to be flipped.'))
            }
            if (assay == 'chipseq'){
                message('ChIP-Seq flip/unflip')
                i <- which(private$table$strand == '-')
                private$table$bin[i] <- (self$get_params()$bin_count + 1) - 
                                                        private$table$bin[i]
                private$table$bin <- as.integer(private$table$bin)
                private$df <- NULL
            } else if (assay == 'rnaseq'){
                message('RNA-Seq flip/unflip')
                i <- which(private$table$strand == '-')
                #col_nuc
                private$table$nuc[i] <- (private$table$regionsize[i] + 1) - 
                                                private$table$nuc[i]
                private$table$nuc <- as.integer(private$table$nuc)
                #col_nuctot
                private$table$nuctot[i] <- (private$table$regionsize[i] + 1) - 
                                        private$table$nuctot[i] + 
                                        private$table$regionstartnuc[i] * 2 - 2
                private$table$nuctot <- as.integer(private$table$nuctot)
                #col_bin
                if(!is.null(private$ph$get("bin_count"))){
                    private$table$bin[i] <- (self$get_params()$bin_count + 1) - 
                                                        private$table$bin[i]
                    private$table$bin <- as.integer(private$table$bin)
                } else {
                    private$table$bin = private$table$nuc
                }
                private$df <- NULL
            }
        },
        get_bam_names = function(filenames) {
            if (all(filenames %in% colnames(private$ph$get("design"))[-1])) {
                filenames
            } else {
                stopifnot(private$check_bam_files(filenames))
                vapply(filenames,
                    private$bam_handler$get_bam_name,
                    character(1))
            }
        },
        check_bam_files = function(bam_files) {
            all(vapply(bam_files,
            function(x) {
                !is.null((private$bam_handler$get_bam_name(x)))
            },
            logical(1)))
        },
        data_frame_avoid_gaps = function(input_df, bam_name, gaps_threshold, flip_regions, bin_count) {
            #bootstrap not executed at this point. Don't work on design !
            
            #how_namy_by_exon_by_design
            dfdt <- data.table::copy(private$df)
            work_df = data.table::copy(input_df)
            nb_nuc_removed <- dfdt[value <= gaps_threshold 
                                    & bam == bam_name, length(value),
                                by=c('exon', 'region')]
            
            #assignment of new exonsize
            for (i in 1:length(nb_nuc_removed$V1)){
                #selected = lines of the ith region and exon of nb_nuc_removed
                selected <- which(
                    work_df$region == nb_nuc_removed$region[i] &
                    work_df$exon == nb_nuc_removed$exon[i])
                #retrieve the exonsize value of the ith region and exon
                original_exonsize <- unique(work_df$exonsize[selected])
                #replace former exonsixe
                new_exonsize <- original_exonsize-nb_nuc_removed$V1[i]
                work_df$exonsize[selected] <- new_exonsize
            }
            
            nb_nuc_removed_by_gene <- dfdt[value <= gaps_threshold 
                                    & bam == bam_name, length(value),
                                by=c('region')]
            #assignment of new region/genesize
            for (i in 1:length(unique(nb_nuc_removed_by_gene$region))){
                #selected = lines of the ith region of nb_nuc_removed
                selected <- which(
                    work_df$region == nb_nuc_removed_by_gene$region[i])
                #retrieve the regionsize value of the ith region and exon
                original_regionsize <- unique(work_df$regionsize[selected])
                #replace former regionsize
                new_regionsize <- (original_regionsize
                                    - nb_nuc_removed_by_gene$V1[i])
                work_df$regionsize[selected] <- new_regionsize
            }
            
            ### removal of zero values
            ## stop if all bam haven't the same amount of lines in table
            stopifnot(length(unique(private$table[, .N, by=bam]$N)) == 1)
            bam_line_count <- private$table[bam == bam_name, .N]
            #lines_to_remove for bam_name
            lines_to_remove <- which(work_df$bam == bam_name &
                                            work_df$value <= gaps_threshold)
            # %% provide the idx for the first bam
            lines_to_remove <- (lines_to_remove %% bam_line_count)
            #to avoid 0 if there is a x %% x = 0
            lines_to_remove <- replace(lines_to_remove, 
                                    which(lines_to_remove == 0), 
                                    bam_line_count)
            bam_count <- length(unique(work_df$bam))
            #lines_to_remove for all bam
            lines_to_remove <- unlist(map((0:(bam_count-1)), 
                                ~ lines_to_remove + bam_line_count * .x))
            work_df <- work_df[-lines_to_remove,]

            
            #reinitialization of nuctot before flip in next section to 
            # clear gaps in nuctot number seauence
            work_df$nuctot <- rep(1:length(which(
                            work_df$bam == bam_name)),
                            times = length(unique(work_df$bam)))
            
            #reorder the nuc and nuctot variables
            if(flip_regions){
                flip_by_bam_n_region <- map2(rep(unique(work_df$bam), 
                            each=length(unique(work_df$region))), 
                    rep(unique(work_df$region),
                            times=length(unique(work_df$bam))), 
                    ~which(work_df$bam == .x & work_df$region == .y 
                                & work_df$strand == '-'))
                
                not_empty_idx <- which(map(flip_by_bam_n_region, 
                                                        ~length(.x)) > 0) 
                if (length(not_empty_idx) > 0){
                    map(flip_by_bam_n_region[not_empty_idx],
                                    ~ (work_df$nuc[.x] <- length(.x):1))
                    map(flip_by_bam_n_region[not_empty_idx],
                                    ~ (work_df$nuctot[.x] <- 
                                            max(work_df$nuctot[.x]):
                                                min(work_df$nuctot[.x])))
                }
                
                unflip_by_bam_n_region <- map2(rep(unique(work_df$bam), 
                            each=length(unique(work_df$region))), 
                    rep(unique(work_df$region), 
                            times=length(unique(work_df$bam))), 
                    ~which(work_df$bam == .x & work_df$region == .y 
                                & (work_df$strand == '+' | 
                                    work_df$strand == '*')))
                not_empty_idx <- which(map(unflip_by_bam_n_region, 
                                                        ~length(.x)) > 0) 
                if (length(not_empty_idx) > 0){
                    map(unflip_by_bam_n_region[not_empty_idx],
                                    ~ (work_df$nuc[.x] <- 1:length(.x)))
                    map(unflip_by_bam_n_region[not_empty_idx],
                                    ~ (work_df$nuctot[.x] <- 
                                            min(work_df$nuctot[.x]):
                                                max(work_df$nuctot[.x])))
                }
            } else { #flip_regions == FALSE
                by_bam_n_region <- map2(rep(unique(work_df$bam), 
                            each=length(unique(work_df$region))), 
                    rep(unique(work_df$region), 
                            times=length(unique(work_df$bam))), 
                    ~which(work_df$bam == .x & work_df$region == .y))
                not_empty_idx <- which(map(by_bam_n_region, 
                                                        ~length(.x)) > 0) 
                if (length(not_empty_idx) > 0){
                    map(by_bam_n_region[not_empty_idx], 
                                    ~ (work_df$nuc[.x] <- 1:length(.x)))
                    map(by_bam_n_region[not_empty_idx], 
                                    ~ (work_df$nuctot[.x] <- 
                                            min(work_df$nuctot[.x]):
                                                max(work_df$nuctot[.x])))
                }
            }
            if(!is.null(bin_count)){
                #reinitialization of region/gene_size to be able to rebuild 
                #bin column
                length_by_region_n_bam <- work_df[,length(nuc),
                                                    by=c('region','bam')]$V1
                work_df$regionsize <- rep(length_by_region_n_bam, 
                                            times=length_by_region_n_bam)
                #rebuild the correct bin column
                col_bins <- trunc((work_df$nuc/(work_df$regionsize+1))
                                        *bin_count)+1
                work_df$bin <- as.integer(col_bins)
            }
            
            return(work_df)
        },
        get_design_names = function(design) {
            if(is.null(design)) {
                return(NULL)
            } else {
                return(colnames(design)[-1])
            }
        },
        get_design_number = function(design) {
            if(is.null(design)) {
                return(NULL)
            } else {
                return(ncol(design) - 1)
            }        
        },
        get_bam_in_design = function(design, design_name) {
            private$get_x_in_design(design, design_name, 1)
        },
        get_control_in_design = function(design, design_name) {
            private$get_x_in_design(design, design_name, 2)        
        },
        get_x_in_design = function(design, design_name, value) {
            if(is.null(design)) {
                return(NULL)
            } else {
                return(design$Samples[design[[design_name]] == value])
            }        
        },
        get_bam_by_design = function(design) {
            map(private$get_design_names(design), ~private$get_bam_in_design(design, .x))
        },
        # Function to generate a table for rna-seq data. 
        # This is a workhorse function which assumes that all parameter
        # validations have already been passed.
        #
        # Here we provide an example of the expected results for
        # a design with two groups (d1 and d2) of two and one bam files respectively
        # (b1 and b2 for d1, b3 for d2):
        #  Samples   d1   d2
        #     b1     1    0
        #     b2     1    0
        #     b3     0    1
        #
        # The resulting table holds its values in the following hierarchical order.
        # Values are shown until they start repeating. Values preceded by
        # a <- indicate that they can be directly inferred from the columns 
        # to their left (Example: exon size can be inferred from exon number),
        # and thus have the same cycle.
        #
        # Most columns have self-explanatory names, except three:
        #  - nuc: represents the nucleotide position within the gene/region.
        #  - regionstartnuc: Multiple genes are "linearized" in a virtual nucleotide space.
        #                    Suppose gene 1 has 6 nucleotides, and gene 2 has five.
        #                    Gene 1 occupies the [1,6] span of the linearized space,
        #                    while gene 2 occupies the [7,11] span.
        #                    This column indicates the start position of the gene
        #                    within the linearized space. It would be 1 for gene 1,
        #                    and 7 for gene 2.
        #  - nuctot: The position in the linearized gene space.
        #
        # Note that nuc + regionstartnuc - 1 = nuctot.
        #
        # Here is the table representing the above examples. Note that the table
        # may be post-processed by flip_regions before being returned.
        #
        # design | bam | region | regionsize | strand | regionstartnuc | exon | exonsize | bin | nuc | nuctot | value | 
        #   d1   |  b1 |   r1   |    <- 6    |  <- +  |     <- 1       | r1e1 |   <- 3   |  1  |  1  |   1    |   x   | 
        #   .    |  .  |   .    |    <- 6    |  <- +  |     <- 1       |  .   |   <- 3   |  1  |  2  |   2    |   x   | 
        #   .    |  .  |   .    |    <- 6    |  <- +  |     <- 1       |  .   |   <- 3   |  2  |  3  |   3    |   x   | 
        #   .    |  .  |   .    |    <- 6    |  <- +  |     <- 1       | r1e2 |   <- 3   |  2  |  4  |   4    |   x   | 
        #   .    |  .  |   .    |    <- 6    |  <- +  |     <- 1       |  .   |   <- 3   |  3  |  5  |   5    |   x   | 
        #   .    |  .  |   .    |    <- 6    |  <- +  |     <- 1       |  .   |   <- 3   |  3  |  6  |   6    |   x   |    
        #   .    |  .  |   r2   |    <- 5    |  <- -  |     <- 7       | r2e1 |   <- 5   |  1  |  1  |   7    |   x   | 
        #   .    |  .  |   .    |    <- 5    |  <- -  |     <- 7       |  .   |   <- 5   |  1  |  2  |   8    |   x   | 
        #   .    |  .  |   .    |    <- 5    |  <- -  |     <- 7       |  .   |   <- 5   |  2  |  3  |   9    |   x   | 
        #   .    |  .  |   .    |    <- 5    |  <- -  |     <- 7       |  .   |   <- 5   |  2  |  4  |   10   |   x   | 
        #   .    |  .  |   .    |    <- 5    |  <- -  |     <- 7       |  .   |   <- 5   |  3  |  5  |   11   |   x   | 
        #   .    |  b2 |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   |    
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   |                     
        #   d2   |  b3 |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   | 
        #   .    |  .  |        |            |        |                |      |          |     |     |        |   x   |                     
        produce_rna_table = function(coverages, design, regions, bin_count) {
          
            # Calculate useful variables. Here the word 'gene' = 'region'
            gene_count <- length(regions)
            gene_names <- names(regions)
            exon_lengths <- width(regions)
            gene_lengths <- vapply(exon_lengths, sum, numeric(1))
            
            # The number of rows/nucleotides for a region is the 
            # sum of all exon lengths.
            nuc_per_region = vapply(exon_lengths, sum, numeric(1))
            
            # The number of rows per bam file the total number of nucleotides
            # within all regions.
            row_per_bam = sum(nuc_per_region)
 
            nb_bfile_by_design <- purrr::map_int(private$get_bam_by_design(design), ~length(.x))

            ##### design column #####
            # The number of rows per design is the number of nucleotides per region
            # times the number of bams in that region.
            rows_per_design = nb_bfile_by_design * row_per_bam
            col_design <- rep(private$get_design_names(design), times=rows_per_design) 
            
            ##### bam column #####
            col_bam = rep(unlist(private$get_bam_by_design(design)), each = row_per_bam)
            
            ##### region/gene columns #####
            col_gene <- rep(gene_names, times=gene_lengths)
            col_gene_size <- rep(gene_lengths, times=gene_lengths)
            
            # Strands are presumed to be identical throughout a region, so
            # we only grab the first one.
            strand_per_gene = unlist(lapply( strand(regions), function(x){as.character(x[1])}))
            col_strand <- rep(strand_per_gene, times=gene_lengths)

            gene_length_cum <- c(0, cumsum(gene_lengths)[-length(gene_lengths)])+1
            col_gene_start_nuc <- rep(gene_length_cum, times = gene_lengths)
            
            ##### exon columns #####
            exon_counts <- vapply(regions, length, numeric(1))                
            exon_names <- unlist(map(exon_counts, ~ 1:.x))
            v_exon_lengths = as.vector(unlist(exon_lengths))
            col_exon <- as.vector(rep(exon_names, times=v_exon_lengths))
            
            #useful for flip function            
            col_exon_size <- rep(v_exon_lengths, times=v_exon_lengths)
                        
            ##### nuc/bin columns #####
            col_nuc_tot = 1:row_per_bam
            col_nuc = (col_nuc_tot - col_gene_start_nuc) + 1

            if (!is.null(bin_count)) {
                col_bins <- as.integer(trunc((col_nuc/col_gene_size)*bin_count)+1)
            } else {
                col_bins = col_nuc
            }
                        
            ## col_values
            #NB : lapply(Views...) -> out of limits of view
            grtot <- self$get_regions()
            col_values <- list()
            idx <- 1 #index for col_values list
            idx_sd_loop <- 1 
            for(bam in unlist(private$get_bam_by_design(design))) {
                bam_name = tools::file_path_sans_ext(basename(bam))
                for (i in 1:length(grtot)){
                    gr <- grtot[[i]]
                    sq <- unique(as.character(seqnames(gr)))
                    val <- Views(
                        coverages[[bam_name]][[sq]], 
                        start(gr), 
                        end(gr))
                    col_values[[idx]] <- unlist(lapply(val, as.numeric))
                    idx <- idx + 1
                }
            }
            col_values <- unlist(col_values)
            
            message('produce data table : RNA-Seq')                                        
            return(data.table(region = col_gene,
                              exon = col_exon,
                              bam = col_bam,
                              design = col_design,
                              bin = col_bins,
                              nuc = col_nuc,
                              nuctot = 1:length(col_nuc),
                              exonsize = col_exon_size,
                              regionstartnuc = col_gene_start_nuc,
                              regionsize = col_gene_size,
                              value = col_values,
                              strand = col_strand))
        },
        # Function to generate a table for chip-seq data. 
        # This is a workhorse function which assumes that all parameter
        # validations have already been passed.
        #
        # Here we provide an example of the expected results for two regions
        # (Each with two ranges) and two designs.
        #
        # The resulting table holds its values in the following hierarchical order.
        # Values are shown until they start repeating.
        #
        # Note that "bin" repeats as many times as there are ranges within
        # the given region.
        #
        # Here is the table representing the above examples. Note that the table
        # may be post-processed by flip_regions before being returned.
        #
        # region | design |    bin     | value | strand |
        #   r1   |  d1    |     1      |   x   |   +    |
        #   .    |  .     |     .      |   x   |   +    |
        #   .    |  .     |  bin_count |   x   |   +    |
        #   .    |  .     |     1      |   x   |   -    |   
        #   .    |  .     |     .      |   x   |   -    |
        #   .    |  .     |  bin_count |   x   |   -    |
        #   .    |  d2    |            |   x   |        |
        #   .    |  .     |            |   x   |        |
        #   .    |  .     |            |   x   |        |
        #   .    |  .     |            |   x   |        |
        #   .    |  .     |            |   x   |        |
        #   .    |  .     |            |   x   |        |   
        #   r2   |  d1    |            |   x   |   -    |
        #   .    |  .     |            |   x   |   -    |
        #   .    |  .     |            |   x   |   -    |
        #   .    |  .     |            |   x   |   +    |
        #   .    |  .     |            |   x   |   +    |
        #   .    |  .     |            |   x   |   +    |
        #   .    |  d2    |            |   x   |        |
        #   .    |  .     |            |   x   |        |
        #   .    |  .     |            |   x   |        |
        #   .    |  .     |            |   x   |        |
        #   .    |  .     |            |   x   |        |
        #   .    |  .     |            |   x   |        |
        produce_chip_table = function(coverages, design, regions, bin_count) {
            message('produce data table : ChIP-Seq')
            region_length <- vapply(regions, length, numeric(1))
            col_regions <- names(regions) %>%
                map(~ rep(.x, length(coverages) * bin_count * 
                            region_length[.x])) %>% unlist()

            #col_regions_range <- as.list(regions) %>%
            #                         map(names) %>%
            #                         map(rep, each=bin_count) %>%
            #                         map(rep, times=length(coverages)) %>%
            #                         unlist()
                                     
            col_designs <- map(region_length, ~ rep(names(coverages), 
                                each = bin_count * .x)) %>% unlist
            col_bins <- rep(1:bin_count,
                            length(coverages) * sum(region_length))
            pairs <- expand.grid(colnames(design)[-1], 
                                names(regions), 
                                stringsAsFactors = FALSE)
            col_values <- map2(pairs$Var1, pairs$Var2,
                ~ private$get_subtable(coverages[[.x]], .y, 
                    bin_count)) %>% unlist
            
            #TODO : improve col_strand production
            # Vectorize? Slower than loop with small data set.
            # col_strand = unlist(lapply((lapply(lapply(strand(regions), rep, each=bin_count), rep, length(names(coverages)))), as.vector))
            
            col_strand <- list()
            idx <- 1
            for (region_names in unique(col_regions)){
                col_strand[[idx]] <- rep(rep(
                    as.vector(strand(regions)[[region_names]]),
                    each=bin_count),length(unique(col_designs)))
                idx <- idx + 1
            }
            col_strand <- unlist(col_strand)
            
            return(data.table(region = col_regions,
            #                  region_range = col_regions_range,
                              design = col_designs,
                              bin = col_bins,
                              value = col_values,
                              strand = col_strand))
        },
        get_coverage_names = function(coverages) {
            stopifnot(length(setdiff(names(coverages), c("+", "-", "*")))==0)
            if(!is.null(coverages[['+']])) {
                return(names(coverages[['+']]))
            } else if(!is.null(coverages[['-']])) {
                return(names(coverages[['-']]))
            } else if(!is.null(coverages[['*']])) {
                return(names(coverages[['*']]))
            }
        },
        produce_data_frame_internal = function(input_table, alpha, sample_count,
                            avoid_gaps, gaps_threshold, bam_name,
                            assay, input_design, bin_count) {
            results <- data.table::copy(input_table)
            
            if(assay=='rnaseq' && !is.null(bin_count)) {
                sample_size_columns = quote(.(region, design))
            } else {
                sample_size_columns = quote(.(design))
            }
        
            # Set up bootstrap analysis.
            sample_size <- input_table[bin == 1,][
                                    ,.N, by = eval(sample_size_columns)][
                                    , .(min(N))]
            sample_size <- as.integer(sample_size)
        
            out_cols <- c("value", "qinf", "qsup")
            bootstrap <- function(df) {
                sampling <- matrix(df$value[sample(seq_along(df$value),
                                        sample_size * sample_count,
                                        replace = TRUE)],
                                ncol = sample_size)
                values <- colMeans(sampling)
                res <- quantile(values, c(alpha/2, 1-(alpha/2)))
                res <- c(mean(df$value), res)
                names(res) <- out_cols
                as.list(res)
            }            
        
            skip_bootstrap = FALSE
            if (assay == 'chipseq') {
                message('produce data frame : ChIP-Seq')
        
                bootstrap_cols = quote(.(region, design, bin))
                unique_col = c("region", "design", "bin", "strand")
            } else {
                if(avoid_gaps) {
                    if (!is.null(bam_name)){
                        bam_name = results$bam[1]
                    }
                    results = private$data_frame_avoid_gaps(results, bam_name, gaps_threshold,
                                                            private$ph$get("flip_regions"), 
                                                            private$ph$get("bin_count"))
                }                
            
                if (is.null(bin_count)) {
                    message('produce data frame : RNA-Seq')
                    
                    bootstrap_cols = quote(.(region, design, nuctot))
                    unique_col = c("region", "design", "nuctot")
                    
                    if(all(rowSums(design[,-1, drop=FALSE]) == 1) &
                        all(colSums(design[,-1, drop=FALSE]) == 1)){
                        results$qinf <- results$value
                        results$qsup <- results$value
                        skip_bootstrap = TRUE
                    }
                } else {
                    message('produce data frame : RNA-Seq binned')
        
                    bootstrap_cols = quote(.(design, bin))
                    unique_col = c("design", "bin")
                }
            }
            
            if(!skip_bootstrap) {
                results <- results[,c(out_cols) := bootstrap(.SD), 
                                    by = eval(bootstrap_cols)]
            }
            results <- unique(results, by=unique_col)
            
            results <- as.data.frame(results)
            results$group <- as.factor(paste(results$design, results$region, sep="_"))
            
            return(results)
        },
        get_raw_coverages_internal = function(filenames = NULL) {
            if (is.null(filenames)) {
                private$coverages
            } else {
                stopifnot(is.character(filenames))
                stopifnot(length(filenames) > 0)
                bam_names <- private$get_bam_names(filenames)
                stopifnot(length(bam_names) == length(filenames))
                lapply(private$coverages, function(x) { x[bam_names] })
            }
        },
        get_normalized_coverages_internal = function(filenames = NULL) {
            # Define a function which will normalize coverage for a single
            # BAM file.
            normalize_coverage <- function(filename) {
                count <- private$bam_handler$get_aligned_count(filename)
                weight <- 1 / (count / 1000000)
                for(strand in c("+", "-", "*")) {
                    if(!is.null(coverages[[strand]])) {
                        coverages[[strand]][[filename]] <- coverages[[strand]][[filename]] * weight
                    }
                }
            }
            
            # Get the raw coverages.
            coverages <- private$get_raw_coverages_internal(filenames)
            
            # Calculate normalized coverages in parallel.
            coverage_names <- private$get_coverage_names(coverages)
            private$parallel_job$launch_job(data = coverage_names,
                                            FUN = normalize_coverage)
            for(strand in c("+", "-", "*")) {                                                
                if(!is.null(coverages[[strand]])) {
                    names(coverages[[strand]]) <- coverage_names
                }
            }
            coverages
        },
        produce_strand_table = function(coverages_s, assay, design, regions, noise_removal, bin_count) {
            if(is.null(coverages_s)) {
                return(NULL)
            } else {
                if (assay == 'rnaseq') {
                    return(private$produce_rna_table(coverages_s, design, regions, bin_count))
                } else { # chipseq
                    if (!is.null(noise_removal)) {
                        coverages_s <- private$remove_controls(coverages_s, design)
                    } else {
                        coverages_s <- private$merge_chip(coverages_s, design)
                    }

                    return(private$produce_chip_table(coverages_s, design, regions, bin_count))
                }
            }        
        },        
        validate_design = function(design) {
            private$validate_design_format(design)
            private$validate_design_values(design)
        },
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
        },
        validate_design_values = function(design) {
            # At least one file must be used in the design
            if (sum(rowSums(design[ , -1, drop=FALSE]) > 0) == 0) {
                stop("At least one BAM file must be used in the design.")
            }

            # Check if used bam files exist.
            non_empty_rows = rowSums(design[, -1, drop=FALSE]) > 0
            if (!all(purrr::map_lgl(design$Samples[non_empty_rows], private$check_bam_files))) {
                warning("At least one BAM file does not exist")
            }
        
        },
        validate_alpha = function(alpha) {
            stopifnot(is.numeric(alpha))
            stopifnot(alpha >= 0 & alpha <= 1)
        },
        validate_bin_count = function(bin_count) {
            if (!is.null(bin_count)) {
                if (!is.numeric(bin_count) || bin_count <= 0 || as.integer(bin_count) != bin_count) {
                    stop("bin_count must be NULL or a positive integer")
                }
            }
        },
        validate_sample_count = function(sample_count) {
            stopifnot(is.numeric(sample_count))
            stopifnot(sample_count > 0)
            stopifnot(as.integer(sample_count) == sample_count)
        },
        validate_avoid_gaps = function(avoid_gaps) {
            stopifnot(is.logical(avoid_gaps))
        },
        validate_gaps_threshold = function(gaps_threshold) {
            stopifnot(is.numeric(gaps_threshold))
            stopifnot(gaps_threshold >= 0)
        },
        validate_noise_removal = function(noise_removal) {
            if (!is.null(noise_removal) && !(noise_removal %in% c("NCIS"))) {
                stop('noise_removal must be NA, NULL, or "NCIS".')
            }
        },
        validate_normalization = function(normalization) {
            if (!is.null(normalization) && normalization != "RPM") {
                stop("normalization must be NA, NULL or \"RPM\".")
            }
        },
        validate_flip_regions = function(flip_regions) {
            if (!is.logical(flip_regions)) {
                stop("flip_regions must be a logical.")
            }
        },
        validate_assay = function(assay) {
            # Check parameters validity
            if (!is.character(assay)) {
                stop("verbose must be a character value")
            }
            assayTypeAuthorized <- c('chipseq', 'rnaseq')
            if (!(tolower(assay) %in% assayTypeAuthorized)) {
                stop("assay values must be one of 'chipseq' or 'rnaseq'")
            }
        },
        validate_verbose = function(verbose) {
            if (!is.logical(verbose)) {
                stop("verbose must be a logicial value (TRUE or FALSE)")
            }       
        },
        validate_force_seqlevels = function(force_seqlevels) {
            if (!is.logical(force_seqlevels)) {
                stop(paste("force_seqlevels must be a logicial ",
                            "value (TRUE or FALSE)",sep=""))
            }        
        },
        validate_padding_size = function(padding_size) {
            if (!(is.numeric(padding_size) || is.integer(padding_size)) ||
                padding_size < 0 || as.integer(padding_size) != padding_size) {
                stop("padding_size must be a non-negative integer")
            }       
        },
        validate_cores = function(cores) {
            isBiocParallel = is(cores, "BiocParallelParam")
            isInteger = ((is.numeric(cores) || is.integer(cores)) &&
                            cores > 0 && as.integer(cores) == cores)
            if (!isBiocParallel && !isInteger) {
                stop("cores must be a positive integer or a BiocParallelParam instance.")
            }        
        },
        validate_bam_files = function(bam_files) {
            private$validate_bam_files_format(bam_files)
            private$validate_bam_files_values(bam_files)
        },
        validate_bam_files_format = function(bam_files) {
            if (!is.vector(bam_files, "character")) {
                stop("bam_files must be a vector of BAM filenames.")
            }      
        },
        validate_bam_files_values = function(bam_files) {
            if (!all(sapply(bam_files, file.exists))) {
                stop("At least one BAM file does not exist.")
            }        
        },        
        validate_combination = function(params) {
            if(params$assay=="chipseq" && is.null(params$bin_count)) {
                stop("bin_count cannot be NULL in chipseq assays.")
            }
        },
        get_complete_design = function(bam_files) {
            bam_files = private$name_from_path(bam_files)
            
            # Concatenate the bam names and the identity matrix, then
            # rename all columns but the first.
            design <- cbind(data.frame(bam_files = bam_files), diag(length(bam_files)))
            colnames(design)[-1] = names(bam_files)
            
            design
        },
        name_from_path = function(file_paths) {
            alt_names = tools::file_path_sans_ext(basename(file_paths))
            if(is.null(names(file_paths))) {
                names(file_paths) = alt_names
            } else {
                ifelse(names(file_paths)=="", alt_names, names(file_paths))
            }
            return(file_paths)
        },
        import_regions_from_disk = function(file_names) {
            file_names = private$name_from_path(file_names)
            import_file <- function(region) {
                ext <- tolower(tools::file_ext(region))
                if (ext == "narrowpeak") {
                    extraCols <- c(signalValue = "numeric",
                                    pValue = "numeric", qValue = "numeric",
                                    peak = "integer")
                    rtracklayer::import(region, format = "BED",
                                        extraCols = extraCols)
                } else if (ext == "broadpeak") {
                    extraCols <- c(signalValue = "numeric",
                                    pValue = "numeric", qValue = "numeric")
                    rtracklayer::import(region, format = "BED",
                                        extraCols = extraCols)
                } else if (ext == "gtf" | ext == "gff") {
                    split(rtracklayer::import(region), 
                            rtracklayer::import(region)$gene_id)
                } else {
                    rtracklayer::import(region)
                }
            }
            regions <- private$parallel_job$launch_job(data = file_names,
                                                    FUN = import_file)
            names(regions) <- names(file_names)
            return(regions)
        },
        deep_clone = function(name, value) {
            # With x$clone(deep=TRUE) is called, the deep_clone gets invoked once for
            # each field, with the name and value.
            if (name == "ph") {
                 # `a` is an environment, so use this quick way of copying
                 value$clone()
            } else {
                # For all other fields, just return the value
                value
            }
        },
        clean_design = function(design, bam_files) {
            # If no design is provided, use the default one.
            if(is.null(design)) {
                design = private$get_complete_design(bam_files)
            }
            
            # NA will be overwritten with the previous design later on.
            if(!is.data.frame(design) && is.na(design)) {
                return(NA)
            }
            
            # Make sure the design is in the correct format (data-frame
            # with at least 2 columns) before we try improving it.
            private$validate_design_format(design)
            
            # Standardize names used in first column to match those
            # used in bam_files.

            design_bams = as.character(design[,1])
            design[,1] = design_bams
            if(all(design_bams %in% names(bam_files))) {
                return(design)
            } else {
                inferred_names = names(private$name_from_path(design_bams))
                if(all(inferred_names %in% names(bam_files))) {
                    design[,1] = inferred_names
                    warning("Modifying first column of design to match the names of bam_files rather than their file name.")
                    return(design)
                } else {
                    stop("Design contains samples absent from the list of bam files provided on initialization.")
                }
            }
        }
    )
)

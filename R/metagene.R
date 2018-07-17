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
                                paired_end_strand_mode=2,
                                region_mode="auto", region_metadata=NULL, ...) {

            # Validate the format of bam_files, since it is used to preprocess certain
            # parameters before initialization.
            validate_bam_files_format(bam_files)
            
            if(region_mode=="auto") {
                region_mode = ifelse(assay=='rnaseq', "stitch", "separate")
            }
            
            # Define default design.
            default_design = private$get_complete_design(bam_files)
            design_metadata = data.frame(design=default_design[,1])
            
            # Initialize parameter handler.
            private$ph <- parameter_manager$new(
                param_values=list(
                    design=default_design,
                    design_metadata=design_metadata,
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
                    alpha=0.05,
                    sample_count=1000,
                    region_mode=region_mode,
                    split_by="region_name",
                    region_filter=TRUE,
                    design_filter=TRUE,
                    resampling_strategy="bin"),
                param_validations=list(
                    design=private$validate_design,
                    bam_files=validate_bam_files,
                    padding_size=validate_padding_size,
                    verbose=validate_verbose,
                    force_seqlevels=validate_force_seqlevels,
                    assay=validate_assay,
                    normalization=validate_normalization,
                    noise_removal=validate_noise_removal,
                    avoid_gaps=validate_avoid_gaps,
                    gaps_threshold=validate_gaps_threshold,
                    bin_count=validate_bin_count,
                    alpha=validate_alpha,
                    sample_count=validate_sample_count,
                    region_mode=validate_region_mode),
                overall_validation=validate_combination)
                
            # Update any other parameter passed as a ... argument.
            # Since the parameter manager is locked, any non-existant
            # parameter will cause an error.
            private$ph$update_params(...)
                
            # Prepare objects for parralel processing.
            validate_cores(cores)
            private$parallel_job <- Parallel_Job$new(cores)
            
            # Prepare bam files
            bm = private$start_bm("Prepare bam files")
            private$bam_handler <- Bam_Handler$new(private$ph$get("bam_files") , cores = cores,
                                        paired_end = paired_end,
                                        strand_specific=strand_specific,
                                        paired_end_strand_mode=paired_end_strand_mode)
            private$stop_bm(bm)

            # Prepare regions
            private$print_verbose("Prepare regions...")
            private$regions <- private$prepare_regions(regions, private$ph$get("region_mode"), region_metadata)

            # Parse bam files
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
        get_regions = function() {
            return(private$regions)
        },
        get_split_regions = function() {
            return(split_regions(private$regions,
                                 private$region_metadata, 
                                 private$ph$get("split_by")))
        },
        get_matrices = function() {
            return(private$split_coverage_cache)
        },
        get_data_frame = function(region_names = NULL, design_names = NULL) {
            if(!is.null(private$ci_meta_df)) {
                results = private$ci_meta_df
            } else if(!is.null(private$ci_df)) {
                results = private$ci_df
            } else {
                return(NULL)
            }
            
            region_subset = TRUE
            if (!is.null(region_names)) {
                stopifnot(is.character(region_names))
                stopifnot(all(region_names %in% unique(private$table$region)))
                region_subset = results$region %in% region_names
            }
            
            design_subset = TRUE
            if (!is.null(design_names)) {
                stopifnot(is.character(design_names))
                stopifnot(all(design_names %in% unique(private$table$design)))
                design_subset = results$region %in% region_names
            }
            return(results[region_subset & design_subset,,drop=F])
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
        produce_data_frame = function(alpha=NA, sample_count=NA, avoid_gaps=NA, 
                                      gaps_threshold=NA, design_metadata=NA) {
            self$add_metadata()
            invisible(self)
        },
        calculate_ci = function(alpha=NA, sample_count=NA, resampling_strategy=NA) {
            # Make sure the previous steps have been completed.
            self$split_coverages_by_regions()
            private$update_params_and_invalidate_caches(alpha, sample_count, resampling_strategy)

            if(is.null(private$ci_df)) {
                bm = private$start_bm("Producing data-frame")
                private$ci_df = calculate_matrices_ci(private$split_coverages,
                                                      private$ph$get('sample_count'), 
                                                      private$ph$get('alpha'),
                                                      private$ph$get('resampling_strategy'),
                                                      private$parallel_job)
                private$stop_bm(bm)
            }
            
            invisible(private$ci_df)
        },
        add_metadata = function(design_metadata=NA) {
            # Make sure the previous steps have been completed.
            self$calculate_ci()
            private$update_params_and_invalidate_caches(design_metadata)
                                                    
            if(is.null(private$ci_meta_df)) {
                filtered_design = private$ph$get("design_metadata")[private$ph$get("design_filter"),, drop=FALSE]
                private$ci_meta_df = add_metadata_to_ci(private$ci_df,
                                                        private$split_metadata_cache,
                                                        filtered_design)
            }
            
            invisible(private$ci_meta_df)
        },
        plot = function(region_names = NULL, design_names = NULL, title = NULL,
                        x_label = NULL, facet_by=NULL, group_by=NULL) {
            # 1. Get the correctly formatted table
            self$add_metadata()

            df <- self$get_data_frame(region_names, design_names)
            
            # 3. Produce the graph
            if (is.null(title)) {
                title <- paste(unique(private$df[["group"]]), collapse=" vs ")
            }
            p <- private$plot_graphic(df = df, title = title, 
                                        x_label = x_label, facet_by=facet_by, group_by=group_by)
            print(p)
            private$graph <- p
            invisible(self)
        },
        export = function(bam_file, region, file) {
            # TODO: Deprecate?
            warning("export is deprecated")
            # region <- private$regions[[region]]
            # param <- Rsamtools::ScanBamParam(which = region)
            # alignments <- GenomicAlignments::readGAlignments(bam_file,
            #                                                     param = param)
            # weight <- 1 - private$bam_handler$get_rpm_coefficient(bam_file)
            # seqlevels(alignments) <- seqlevels(region)
            # # TODO: don't use the weight param of coverage
            # coverage <- GenomicAlignments::coverage(alignments, weight=weight)
            # rtracklayer::export(get_normalized_coverages()[bam_file], file, "BED")
            # invisible(coverage)
        },
        group_coverages = function(design=NA, normalization=NA, noise_removal=NA, design_filter=NA) {
            # Clean up the design so it'll have the expected format.
            design = private$clean_design(design, private$ph$get("bam_files"))

            if(private$ph$have_params_changed(design)) {
                # design has changed.
                # Generate a new design metadata.
                private$ph$set("design_metadata", data.frame(design=design[,1]))
            }
            
            private$update_params_and_invalidate_caches(design, normalization, noise_removal, design_filter)
            
            if(is.null(private$grouped_coverages)) {
                bm <- private$start_bm("Grouping and normalizing coverages")
                design_col_to_keep = rep_len(private$ph$get("design_filter"), ncol(private$ph$get("design")) - 1)
                design_col_to_keep = 1 + which(design_col_to_keep)            
                private$grouped_coverages = lapply(private$get_coverages_internal(),
                                                   group_coverages_s,
                                                   private$ph$get("design")[, c(1, design_col_to_keep)],
                                                   private$ph$get("noise_removal"),
                                                   private$bam_handler)
                private$stop_bm(bm)                                                   
            }

            
            invisible(private$grouped_coverages)
        },
        bin_coverages = function(bin_count=NA, region_filter=NA) {
            # Make sure the previous step has been performed.
            self$group_coverages()
            private$update_params_and_invalidate_caches(bin_count, region_filter)
            
            if(is.null(private$binned_coverages)) {
                bm = private$start_bm("Binning coverages")
                private$binned_coverages = bin_coverages_s(private$grouped_coverages,
                                                           regions=private$regions[private$ph$get("region_filter")],
                                                           bin_count=private$ph$get("bin_count"))
                private$stop_bm(bm)
            }
           
            invisible(private$binned_coverages)
        },
        split_coverages_by_regions = function(split_by=NA) {
            # Make sure the previous step has been performed.
            self$bin_coverages()
            private$update_params_and_invalidate_caches(split_by)
            
            if(is.null(private$split_coverages)) {
                bm = private$start_bm("Splitting coverages by region type")
                split_res = split_matrices(private$binned_coverages,
                                           private$region_metadata[private$ph$get("region_filter"),, drop=FALSE],
                                           private$ph$get('split_by'))
                private$split_coverages = split_res$Matrices
                private$split_metadata_cache = split_res$Metadata
                private$stop_bm(bm)                                         
            }
        },
        produce_metagene = function(...) {
            private$update_params_and_invalidate_caches(...)
            self$plot()
            
            invisible(self)
        }
    ),
    private = list(
        regions = GRangesList(),
        region_metadata = NULL,
        split_regions_cache = NULL,
        split_metadata_cache = NULL,
        design_metadata_cache=NULL,
        coverages = list(),
        grouped_coverages = NULL,
        binned_coverages = NULL,
        split_coverages = NULL,
        df = NULL,
        graph = "",
        bam_handler = "",
        parallel_job = "",
        ph=NULL,
        ci_df=NULL,
        ci_meta_df=NULL,
        
        print_verbose = function(to_print) {
            if (private$ph$get("verbose")) {
                cat(paste0(to_print, "\n"))
            }
        },
        prepare_regions = function(regions, region_mode, region_metadata) {
            if (class(regions) == "character") {
                # Validation specific to regions as a vector
                if (!all(sapply(regions, file.exists))) {
                    stop("regions is a list of files, but some of those files do not exist.")
                }
                regions = private$import_regions_from_disk(regions)
            } else if (class(regions) == "GRanges") {
                if(region_mode=="stitch") {
                    if(length(unique(seqnames(regions))) > 1) {
                        stop(paste0("In stitch region_mode, such as in rnaseq assays, regions should be a ",
                                    "GRangesList of transcripts, or a GRanges ",
                                    " object representing a single transcript. ",
                                    "Here regions spans several seqnames, indicating ",
                                    "it might include many transcripts."))
                    }
                }
                regions <- GRangesList("regions" = regions)                
            } else if (class(regions) == "list") {
                regions <- GRangesList(regions)
            } else if (!is(regions, "GRangesList")) {
                stop(paste0("regions must be either a vector of BED ",
                            "filenames, a GRanges object or a GrangesList object"))            
            }
            
            # If regions do not have names, generate generic names for them.
            if (is.null(names(regions))) {
                names(regions) <- paste0("region_", seq_along(regions))
            }
            
            # In stitch mode, make sure disjoint regions part of the same
            # group do not overlap each other.
            if (region_mode=="stitch"){
                # If some "exons" overlap, then the total size will be smaller than the reduced size.
                total_size = sum(width(regions))
                reduced_size = sum(width(GenomicRanges::reduce(regions)))
                if(!all(total_size==reduced_size)) {
                    stop("In stitch region_mode, no overlap should exist between the individual ",
                         "GRanges making up the elements of the GRangesList")
                }
            }
            # TODO: Check if there is a id column in the mcols of every ranges.
            #    If not, add one by merging seqnames, start and end.

            # Apply padding and sortSeqLevels
            pad_regions = function(x, padding_size) {
                start(x) <- pmax(start(x) - padding_size, 1)
                end(x) <- end(x) + padding_size
                # Clean seqlevels
                x <- sortSeqlevels(x)
                x
            }
            regions = GRangesList(lapply(regions, pad_regions, padding_size = private$ph$get("padding_size")))
            
            # Add a region column to all GRangesList elements.
            regions_with_extra_col = list()
            for(region_name in names(regions)) {
                regions_with_extra_col[[region_name]] = regions[[region_name]]
                mcols(regions_with_extra_col[[region_name]])$region_name = region_name
            }
            regions = GRangesList(regions_with_extra_col)
            
            # In separate mode, simplify regions into a single GRanges object.
            if(region_mode=="separate") {
                regions = unlist(regions, use.names=FALSE)
                private$region_metadata = mcols(regions)
            } else {
                first_or_null = function(x) {
                    if(length(x)>0) {
                        return(mcols(x)[1,, drop=FALSE])
                    } else {
                        return(NULL)
                    }
                }
                private$region_metadata = do.call(rbind, lapply(regions, first_or_null))
            }
            
            return(regions)
        },
        produce_coverages = function() {
            if(private$ph$get("region_mode")=="stitch") {
                regions = BiocGenerics::unlist(private$regions)
            } else {
                regions = private$regions
            }
            
            regions <- GenomicRanges::reduce(regions)
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
        plot_graphic = function(df, title, x_label, facet_by, group_by) {
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
            p <- plot_metagene(df, facet_by=facet_by, group_by=group_by) +
                ylab(y_label) +
                xlab(x_label) +
                ggtitle(title)
            p
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
        get_coverages_internal = function() {
            if (!is.null(private$ph$get("normalization"))) {
                coverages <- private$get_normalized_coverages_internal()    
            } else {
                coverages <- private$get_raw_coverages_internal()
            }
            
            return(coverages)
        },
        validate_design = function(design) {
            validate_design_format(design)
            private$validate_design_values(design)
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
            validate_design_format(design)
            
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
        },
        start_bm = function(msg) {
            private$print_verbose(paste0(msg, "..."))
            return(list(Message=msg, Time=Sys.time(), Memory=pryr::mem_used()))
        },
        stop_bm = function(bm_obj) {
            bm_after_time = Sys.time()
            bm_after_mem = pryr::mem_used()            
            bm_time = difftime(bm_after_time, bm_obj$Time, unit="secs")
            bm_mem = structure(bm_after_mem - bm_obj$Memory, class="bytes")
            
            private$print_verbose(paste0("BENCHMARK-TIME-", bm_obj$Message, ":", bm_time))
            private$print_verbose(paste0("BENCHMARK-MEMORY-", bm_obj$Message, ":", bm_mem))        
        },
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
        },
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
        },
        update_params_and_invalidate_caches = function(...) {
            # This prologue makes it possible to infer parameter names from the
            # name of the variable it is passed in. This allows us to avoid
            # design=design, bin_count=bin_count repetitive code.
            #
            # It cannot be factorized into a function, since in any further call,
            # the argument list will deparse as "list(...)".
            param_names_alt = sapply( substitute(list(...)), deparse)[-1]
            arg_list = list(...)
            if(is.null(names(arg_list))) {
                names(arg_list) = param_names_alt
            } else {
                names(arg_list) = ifelse(names(arg_list)=="", param_names_alt, names(arg_list))
            }
        
            # Associate each parameter witht he step it is used in.
            param_step_map = c(design="group_coverages",
                               normalization="group_coverages", 
                               noise_removal="group_coverages", 
                               design_filter="group_coverages",
                               bin_count="bin_coverages", 
                               region_filter="bin_coverages",
                               split_by="split_coverages",
                               region_filter="split_coverages",
                               alpha="calculate_ci", 
                               sample_count="calculate_ci", 
                               resampling_strategy="calculate_ci",
                               design_metadata="add_metadata")
                               
            # Associate each step with the cache it generates,
            # in reverse order, so we can easily determine which
            # caches to invalidate when a particular step needs to be re-run.
            step_cache_map = c(add_metadata="ci_meta_df",
                               calculate_ci="ci_df",
                               split_coverages="split_coverages",
                               bin_coverages="binned_coverages",
                               group_coverages="grouped_coverages")
                               
            cache_invalidated=FALSE
            
            # Loop over all passed-in parameters.
            for(arg_index in 1:length(arg_list)) {
                arg_name=names(arg_list)[arg_index]
                cat("Analyzing ", arg_name, "\n")
                # Determine if the parameter has changed from its last value.
                if(do.call(private$ph$have_params_changed, arg_list[arg_index])) {
                    cat(arg_name, " has changed.\n")
                    # Determine which step the parameter belongs to.
                    invalidated_step = param_step_map[names(arg_list)[arg_index]]
                    if(!is.na(invalidated_step)) {
                        
                        # Invalidate all caches for the step the parameter belonged to,
                        # as well as all caches for downsteam steps.
                        invalidated_caches = step_cache_map[1:which(names(step_cache_map)==invalidated_step)]
                        cat(paste0(invalidated_caches, collapse=", "), " will be invalidated.\n")
                        for(cache in invalidated_caches) {
                            private[[cache]] = NULL
                            cache_invalidated = TRUE
                        }
                    }
                }
            }
            
            do.call(private$ph$update_params, arg_list)
            
            return(cache_invalidated)
        }        
    )
)

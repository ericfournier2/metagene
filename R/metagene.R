#' A class to manage metagene analysis.
#'
#' This class will allow to load, convert and normalize alignments and regions
#' files/data. Once the data is ready, the user can then choose to produce
#' metagene plots on the data or some subset of it.
#'
#' @section Constructor:
#' \describe{
#'    \item{}{\code{mg <- metagene2$new(regions, bam_files, padding_size = 0,
#'                                     cores = SerialParam(), verbose = FALSE,
#'                                     force_seqlevels = FALSE, paired_end = FALSE,
#'                                     assay = 'chipseq', strand_specific=FALSE,
#'                                     paired_end_strand_mode=2,
#'                                     region_mode="auto", region_metadata=NULL, ...))}}
#'    \item{regions}{A description of all regions over which metagenes will be calculated.
#'
#'                   When \code{region_mode} is "separate", those should be provided using
#'                   a \code{GRanges} object representing all individual, contiguous regions
#'                   to be examined.
#'
#'                   When \code{region_mode} is "stitch", those should be provided using a 
#'                   \code{GRangesList} object where each individual \code{GRanges} element
#'                   represents a set of regions to be stitched together.
#'
#'                   As a convenience, in "separate" mode, \code{metagene} will convert any
#'                   passed in \code{GRangesList} into an unlisted \code{GRanges} with an 
#'                   additional region_name metadata column containing the name of the
#'                   \code{GRangesList} element it was extracted from.
#'
#'                   Also as a convenience, regions can also be a \code{character} 
#'                   \code{vector} of filenames, which are then imported into a GRangesList.
#'                   Supported file formats are BED, narrowPeak, broadPeak, gff and gtf.}
#'    \item{bam_files}{A \code{vector} of BAM filenames. The BAM files must be
#'                    indexed. i.e.: if a file is named file.bam, there must
#'                    be a file named file.bam.bai or file.bai in the same 
#'                    directory. If bam_files is a named vector, then the provided names
#'                    can be used downstream to refer to those bam files. If no
#'                    names are provided, \code{metagene} will try to infer appropriate ones.}
#'    \item{assay}{\code{'chipseq'}, \code{'rnaseq'} or NULL. If non-NULL, metagene will
#'                 set other parameters, such as region_mode and strand_specific, to logical
#'                 values for the given assay. Default: \code{'chipseq'}}
#'    \item{region_mode}{Set the way the \code{regions} parameter is interpreted. Can be 
#'                       \code{'separate'}, \code{'stitch'} or \code{'auto'}. In separate mode,
#'                       \code{regions} is expected to be a GRanges defining individual, contiguous regions. In
#'                       \code{'stitch'} mode, \code{regions} is expected to be a GRangesList 
#'                       where each \code{GRanges} element represents a set of regions to be 
#'                       stitched together and treated as a single logical region. If \code{'auto'}
#'                       then a logical value is inferred from the \code{assay} parameter.
#'                       Default: \code{'auto'}}
#'    \item{region_metadata}{A data-frame of metadata to be associated with the elements of regions
#'                           it must contain has many rows as there are elements in regions. If 
#'                           \code{region_metadata} is NULL but \code{regions} has an mcols element,
#'                           then it is used.}
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
#'    \item{paired_end}{Set this to \code{TRUE} if the provided bam files
#'                      describe paired-end reads. If \code{FALSE}, single-ended
#'                      data are expected. Default: \code{FALSE}}
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
#'    \item{...}{Additional parameters for the metagene analysis. See \code{produce_metagene}
#'               for a list of possible parameters.}
#' }
#'
#'    \code{metagene2$new} returns a \code{metagene} object that contains the
#'        coverages for every BAM files in the regions from the \code{regions}
#'        param.
#'
#' @return
#' \code{metagene2$new} returns a \code{metagene} object which contains the
#' normalized coverage values for every regions in all specified BAM files.
#'
#' @section Methods:
#' \describe{
#'    \item{}{\code{mg$group_coverages(design=NA, normalization=NA, 
#'            design_filter=NA)}}
#'    \item{Description}{This method normalizes genome-wide coverages, then groups
#'            them according to the specified design groups. It returns
#'            a list of possible read orientations (+, -, *), each element
#'            of which is either NULL (depending on the value of the 
#'            strand_specific parameter) or a list of possible design groups.
#'            In turn, the lists of design groups contain lists of \code{Rle}
#'            objects representing coverage over a specific chromosome or sequence.}
#'    \item{design}{A \code{data.frame} that describe to experiment to plot.
#'            see \code{plot} function for more details. \code{NA} can 
#'            be used keep previous design value. Default: \code{NA}.}
#'    \item{normalization}{The algorithm to use to normalize samples. Possible
#'                        values are \code{NULL}, "RPM" and "NCIS". See
#'                        Liand and Keles 2012 for the NCIS algorithm.
#'                        Default: NULL.}
#'    \item{design_filter}{}
#' }
#' \describe{
#'    \item{}{\code{mg$bin_coverages(bin_count=NA, region_filter=NA)}}
#'    \item{}{This method summarizes the coverage over regions of interests
#'            into a specified number of bins. It returns a list where each element
#'            represents a design group, and contains a matrix of binned coverages
#'            where each row represents a region, and each column represents a bin.}
#'    \item{bin_count}{The number of bin to create. \code{NA} can be used to
#'                    keep previous bin_count value. A bin_count value of 100
#'                    will be used if no value is specified. Default:
#'                    \code{NA}.}
#'    \item{region_filter}{}
#' }
#' \describe{
#'    \item{}{\code{mg$split_coverages_by_regions(split_by=NA)}}
#'    \item{}{This methods splits the matrices generated by mg$bin_coverages
#'            into groups of regions where the values of the metadata columns
#'            specified by \code{split_by} are homogeneous. It returns A list
#'            where each element represents a design group: each of those
#'            element is in turn a list representing groups of regions for which
#'            all metadata values specified by "split_by" are equal. The leaf elements
#'            of this list hierarchy are coverage matrices where each row represents a
#'            region, and each column represents a bin.}
#'    \item{split_by}{A vector of column names from the region_metadata
#'                    parameter, as specified on metagene initialization. The 
#'                    selected columns must allow conversion into a factor.
#'                    By default, this is set to region_name, a metadata column
#'                    which is automatically generated by metagene.}
#' }
#' \describe{
#'    \item{}{\code{mg$add_metadata(design_metadata=NA)}}
#'    \item{design_metadata}{A data-frame containing metadata for the design groups.
#'                           It must contain as many rows as there are design groups,
#'                           and must contain at least one column named 'design'
#'                           which is used to match the rows to design groups.}
#' }
#' \describe{
#'    \item{}{\code{mg$plot(region_names = NULL, design_names = NULL,
#'                title = NA, x_label = NA, facet_by=NA, group_by=NA)}}
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
#'    \item{}{\code{mg$calculate_ci(alpha = NA, sample_count = NA, resampling_strategy=NA)}}
#'    \item{alpha}{The alpha level of the confidence interval estimate.
#'                \code{1 - alpha / 2} and \code{alpha / 2} will be used.
#'                Default: 0.05.}
#'    \item{sample_count}{The number of draws to perform in the bootstrap
#'                        calculations. Default: 1000.}
#'    \item{resampling_strategy}{The resampling strategy to be sued when performing the
#'                               bootstrap analysis, which can be either \code{'profile'}
#'                               or \code{'bin'}. In \code{'profile'} mode, whole profiles
#'                               across all bins are resampled. In \code{'bin'} mode,
#'                               each bin is resampled individually and independantly from
#'                               all others. }
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
#' mg <- metagene2$new(regions = region, bam_files = bam_file)
#' \dontrun{
#'    df <- metagene2$plot()
#' }
#'
#' @importFrom R6 R6Class
#' @importFrom data.table data.table
#' @export
#' @format A metagene experiment manager

metagene2 <- R6Class("metagene",
    public = list(
    # Methods
        initialize = function(regions, bam_files, padding_size = 0,
                                cores = SerialParam(), verbose = FALSE,
                                force_seqlevels = FALSE, paired_end = FALSE,
                                assay = 'chipseq', strand_specific=FALSE,
                                paired_end_strand_mode=2,
                                region_mode="auto", region_metadata=NULL, 
                                extend_reads=0, ...) {

            # Validate the format of bam_files, since it is used to preprocess certain
            # parameters before initialization.
            validate_bam_files_format(bam_files)
            
            if(region_mode=="auto") {
                region_mode = ifelse(assay=='rnaseq', "stitch", "separate")
            }
            
            # Define default design.
            default_design = private$get_complete_design(bam_files)
            design_metadata = data.frame(design=as.character(default_design[,1]), stringsAsFactors=FALSE)
            
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
                    avoid_gaps=FALSE,
                    gaps_threshold=0,
                    bin_count=100,
                    alpha=0.05,
                    sample_count=1000,
                    region_mode=region_mode,
                    split_by="region_name",
                    region_filter=TRUE,
                    design_filter=TRUE,
                    resampling_strategy="bin",
                    facet_by=NULL,
                    group_by=NULL,
                    title=NULL,
                    x_label=NULL,
                    extend_reads=extend_reads),
                param_validations=list(
                    design=private$validate_design,
                    bam_files=validate_bam_files,
                    padding_size=validate_padding_size,
                    verbose=validate_verbose,
                    force_seqlevels=validate_force_seqlevels,
                    assay=validate_assay,
                    normalization=validate_normalization,
                    avoid_gaps=validate_avoid_gaps,
                    gaps_threshold=validate_gaps_threshold,
                    bin_count=validate_bin_count,
                    alpha=validate_alpha,
                    sample_count=validate_sample_count,
                    region_mode=validate_region_mode,
                    extend_reads=validate_extend_reads),
                overall_validation=validate_combination)
                
            # Update any other parameter passed as a ... argument.
            # Since the parameter manager is locked, any non-existant
            # parameter will cause an error.
            private$ph$update_params(...)
                
            # Prepare objects for parralel processing.
            self$set_cores(cores)
            
            # Prepare bam files
            bm = private$start_bm("Prepare bam files")
            private$bam_handler <- Bam_Handler$new(private$ph$get("bam_files") , cores = cores,
                                        paired_end = paired_end,
                                        strand_specific=strand_specific,
                                        paired_end_strand_mode=paired_end_strand_mode,
                                        extend_reads=private$ph$get("extend_reads"))
            private$stop_bm(bm)

            # Prepare regions
            private$print_verbose("Prepare regions...")
            private$regions <- private$prepare_regions(regions, private$ph$get("region_mode"), region_metadata)

            # Parse bam files
            bm = private$start_bm("Producing coverage")
            private$coverages <- private$produce_coverages()
            private$stop_bm(bm)            
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
        get_design_group_names = function() {
            private$get_design_names_internal(private$ph$get('design'))
        },        
        get_regions = function() {
            return(private$regions)
        },
        get_regions_metadata = function() {
            return(private$region_metadata)
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
                region_subset = results$region %in% region_names
            }
            
            design_subset = TRUE
            if (!is.null(design_names)) {
                design_subset = results$design %in% design_names
            }
            return(results[region_subset & design_subset,,drop=F])
        },
        get_plot = function() {
            private$graph
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
        set_cores = function(cores) {
            validate_cores(cores)
            if(is.null(private$parallel_job)) {
                private$parallel_job <- Parallel_Job$new(cores)
            } else {
                private$parallel_job$set_core_count(cores)
            }
        },
        produce_data_frame = function(alpha=NA, sample_count=NA, avoid_gaps=NA, 
                                      gaps_threshold=NA, design_metadata=NA) {
            self$add_metadata()
            invisible(self)
        },
        group_coverages = function(design=NA, normalization=NA, design_filter=NA) {
            # Clean up the design so it'll have the expected format.
            design = private$clean_design(design, private$ph$get("bam_files"))

            if(private$ph$have_params_changed(design)) {
                # design has changed.
                # Generate a new design metadata.
                private$ph$set("design_metadata", data.frame(design=design[,1]))
            }
            
            private$update_params_and_invalidate_caches(design, normalization, design_filter)
            
            if(is.null(private$grouped_coverages)) {
                bm <- private$start_bm("Grouping and normalizing coverages")
                
                # Identify design subset.
                design_col_to_keep = rep_len(private$ph$get("design_filter"), ncol(private$ph$get("design")) - 1)
                design_col_to_keep = 1 + which(design_col_to_keep)
                design_subset = private$ph$get("design")[, c(1, design_col_to_keep)]
                
                # Determine how data is to be merged
                if(is.null(private$ph$get("normalization"))) {
                    merge_operation = "+"
                } else if(private$ph$get("normalization")=="RPM") {
                    merge_operation = "mean"
                } else if(private$ph$get("normalization")=="NCIS") {
                    merge_operation = "NCIS"
                } else {
                    stop("Unsupported normalization value.")
                }
                
                private$grouped_coverages = lapply(private$get_coverages_internal(),
                                                   group_coverages_s,
                                                   design=design_subset,
                                                   bam_handler=private$bam_handler,
                                                   merge_operation=merge_operation)
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
                                                           regions=private$select_regions(private$ph$get("region_filter")),
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
                                           private$select_region_metadata(private$ph$get("region_filter")),
                                           private$ph$get('split_by'))
                private$split_coverages = split_res$Matrices
                private$split_metadata_cache = split_res$Metadata
                private$stop_bm(bm)                                         
            }
            
            invisible(private$split_coverages)
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
        plot = function(region_names = NULL, design_names = NULL, title = NA,
                        x_label = NA, facet_by=NA, group_by=NA) {
            # 1. Get the correctly formatted table
            self$add_metadata()

            private$update_params_and_invalidate_caches(title, x_label, facet_by, group_by)
            
            plot_df <- self$get_data_frame(region_names, design_names)
            
            # 3. Produce the graph
            if (is.null(title)) {
                title <- paste(unique(plot_df[["group"]]), collapse=" vs ")
            }
            private$graph <- private$plot_graphic(df = plot_df, 
                                        title = private$ph$get("title"), 
                                        x_label = private$ph$get("x_label"),
                                        facet_by=private$ph$get("facet_by"),
                                        group_by=private$ph$get("group_by"))
            private$graph
        },
        produce_metagene = function(...) {
            private$update_params_and_invalidate_caches(...)
            self$plot()
        },
        plot_single_region = function(region, facet_by=NA, group_by=NA,
                                      no_binning=FALSE) {
            # Clone the mg object
            single_mg = self$clone(deep=TRUE)
            
            # Select the one region to plot, and make sure it ends up as a single GRange object.
            single_region = private$select_regions(region)
            stopifnot(length(single_region)==1)
            if(class(self$get_regions()) == "GRangesList") {
                single_region = single_region[[1]]
            }
            
            # If we are skipping binning, make the number of bins equal
            # the length in nucleotide.
            if(no_binning) {
                bin_count = sum(width(single_region))
            } else {
                bin_count = private$ph$get("bin_count")
            }
            
            # Re-bin with the new bin_count and the new filter.
            single_mg$bin_coverages(bin_count=bin_count, region_filter=region)
            
            # With a single region, splitting cannot be applied, so
            # we'll pass in a single row-name that we know to be valid.
            single_mg$split_coverages_by_regions(split_by="region_name")
            
            # Produce the new base single plot.
            out_plot = single_mg$plot(facet_by=facet_by, group_by=group_by)
            
            # Adjust the plot if binning was skipped.
            if(no_binning) {
                # New x label.
                out_plot <- out_plot + labs(x="Distance in nucleotides") 
                
                # If dealing with a stitched region, display "exon" boundaries as.
                # vertical lines.
                if(length(single_region) > 1) {
                    cumulative_width = 0
                    for(i in width(single_region)) {
                        cumulative_width = cumulative_width + i
                        out_plot <- out_plot + geom_vline(xintercept=cumulative_width)
                    }
                }
            }

            # Return the new plot.
            out_plot
        },
        replace_region_metadata = function(region_metadata) {
            # Validate that the old and new metadata have the same number of rows.
            if(nrow(region_metadata)!=nrow(private$region_metadata)) {
                stop("region_metadata must have one row per region.")
            }
            
            # Make sure the region_name column is still present.
            if(is.null(region_metadata$region_name)) {
                warning("region_name is missing from the new metadata. Recreating it.")
                region_metadata$region_name = private$region_metadata$region_name
            }
            
            # Make sure that the split_by columns are all present and did not change.
            # If they did, invalidate everything after split_by
            for(split_column in private$ph$get("split_by")) {
                if(is.null(region_metadata[[split_column]]) ||
                   !all(region_metadata[[split_column]]==private$region_metadata[[split_column]])) {
                    # The split_by columns were removed or changed. Restore the original parameter value.
                    warning("Replace region_metadata with metadata which would result in a different ",
                            "region split. All caches at the 'split_regions' step will be invalidated, ",
                            "and split_by will be reset to its default value.")
                    private$update_params_and_invalidate_caches(split_by="region_name")
                }
            }
            
            # Everything checks out, replace the metadata.
            private$region_metadata = region_metadata
        }
    ),
    private = list(
        # Region information. Both are kept separate from the parameter
        # handler since both can be very large, and making comparisons
        # to see if they've changed would be onerous.
        regions = GRangesList(),
        region_metadata = NULL,
        
        # Caches of intermediary step results.
        coverages = list(),
        grouped_coverages = NULL,
        binned_coverages = NULL,
        split_coverages = NULL,
        ci_df=NULL,
        ci_meta_df=NULL,
        graph = NULL,

        # Internal caches.
        split_metadata_cache = NULL,

        # Objects for handling bams, parameters and parallel jobs.
        bam_handler = "",
        parallel_job = NULL,
        ph=NULL,
        
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
            }

            if (class(regions) == "GRanges") {
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
            }

            # Build metadata from the mcols of the given regions.
            if(region_mode=="separate") {
                mcol_metadata = mcols(regions)
            } else {
                first_or_null = function(x) {
                    if(length(x)>0) {
                        return(mcols(x)[1,, drop=FALSE])
                    } else {
                        return(NULL)
                    }
                }
                mcol_metadata = do.call(rbind, lapply(regions, first_or_null))
                rownames(mcol_metadata) = names(regions)
            }

            # Merge the passed metadata object with the mcol metadata.
            if(is.null(region_metadata)) {
                private$region_metadata = mcol_metadata
            } else {
                stopifnot(nrow(region_metadata)==length(regions))
                non_duplicate_columns = setdiff(colnames(mcol_metadata), colnames(region_metadata))
                if(length(non_duplicate_columns) > 0) {
                    private$region_metadata = cbind(region_metadata, mcol_metadata[,non_duplicate_columns, drop=F])
                }
            }
            
            if(!is.null(names(regions)) && is.null(rownames(private$region_metadata))) {
                rownames(private$region_metadata) = names(regions)
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

            names(res) <- bam_files
            
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
        get_design_names_internal = function(design) {
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
            map(private$get_design_names_internal(design), ~private$get_bam_in_design(design, .x))
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
                stopifnot(all(filenames %in% names(private$coverages)))
                return(private$coverages[filenames])
            }
        },
        get_normalized_coverages_internal = function(filenames = NULL) {
            # Define a function which will normalize coverage for a single
            # BAM file.
            normalize_coverage <- function(work_item) {
                weight <- 1 / (work_item$Count / 1000000)
                return(list(Strand=work_item$Strand, 
                            BamFile=work_item$BamFile, 
                            Coverage=work_item$Coverage * weight))
            }
            
            # Get the raw coverages.
            coverages <- private$get_raw_coverages_internal(filenames)
            
            # Serialize the workload
            work_items = list()
            for(strand in c("+", "-", "*")) {
                if(!is.null(coverages[[strand]])) {
                    for(bam_file in names(coverages[[strand]])) {
                        work_items[[length(work_items) + 1]] = list(Coverage=coverages[[strand]][[bam_file]],
                                                                  Strand=strand,
                                                                  BamFile=bam_file,
                                                                  Count=private$bam_handler$get_aligned_count(bam_file))
                    }
                }
            }
            
            serialized_coverages <- private$parallel_job$launch_job(data = work_items,
                                                                    FUN = normalize_coverage)
                                                                    
            # Deserialize the coverages.
            norm_coverages = list("+"=NULL, "-"=NULL, "*"=NULL)
            for(i in serialized_coverages) {
                if(is.null(norm_coverages[[i$Strand]])) {
                    norm_coverages[[i$Strand]] = list()
                }
                
                norm_coverages[[i$Strand]][[i$BamFile]] = i$Coverage
            }

            norm_coverages
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
            if (!all(file.exists(as.character(design[non_empty_rows,1])))) {
                warning("At least one BAM file does not exist")
            }
        },        
        get_complete_design = function(bam_files) {
            bam_files = private$name_from_path(bam_files)
            
            # Concatenate the bam names and the identity matrix, then
            # rename all columns but the first.
            design <- cbind(data.frame(bam_files = bam_files, stringsAsFactors=FALSE), diag(length(bam_files)))
            colnames(design)[-1] = names(bam_files)
            
            design
        },
        name_from_path = function(file_paths) {
            alt_names = tools::file_path_sans_ext(basename(file_paths))
            if(is.null(names(file_paths))) {
                names(file_paths) = alt_names
            } else {
                names(file_paths) = ifelse(names(file_paths)=="", alt_names, names(file_paths))
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
                    gxf_regions = rtracklayer::import(region)
                    if("gene_id" %in% colnames(mcols(gxf_regions))) {
                        return(split(gxf_regions, gxf_regions$gene_id))
                    } else {
                        return(gxf_regions)
                    }
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
            if(!all(as.character(design[,1]) %in% bam_files)) {
                stop("Design contains samples absent from the list of bam files provided on initialization.")
            }

            return(design)
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
            if(length(arg_list) > 0) {
                for(arg_index in 1:length(arg_list)) {
                    arg_name=names(arg_list)[arg_index]
                    
                    # Determine if the parameter has changed from its last value.
                    if(do.call(private$ph$have_params_changed, arg_list[arg_index])) {
                        private$print_verbose(paste0(arg_name, " has changed.\n"))
                        
                        # Determine which step the parameter belongs to.
                        invalidated_step = param_step_map[names(arg_list)[arg_index]]
                        if(!is.na(invalidated_step)) {
                            # Invalidate all caches for the step the parameter belonged to,
                            # as well as all caches for downsteam steps.
                            invalidated_caches = step_cache_map[1:which(names(step_cache_map)==invalidated_step)]
                            private$print_verbose(paste0(paste0(invalidated_caches, collapse=", "), " will be invalidated.\n"))
                            for(cache in invalidated_caches) {
                                private[[cache]] = NULL
                                cache_invalidated = TRUE
                            }
                        }
                    }
                }
            
                do.call(private$ph$update_params, arg_list)
            }
            
            return(cache_invalidated)
        },
        select_region_indices = function(selector) {
            if("quosure" %in% class(selector)) {
                # Using dplyr for this because I'm not comfortable enough with
                # quosures.
                selected_indices = as.data.frame(private$region_metadata) %>%
                                       dplyr::mutate(METAGENE_IDX=1:n()) %>%
                                       dplyr::filter(!! selector) %>%
                                       dplyr::pull(METAGENE_IDX)
            } else if(class(selector)=="character") {
                if(!is.null(names(self$get_regions()))) {
                    selected_indices = selector
                } else {
                    selected_indices = self$get_regions()$region_name %in% selector
                }
            } else if(is.numeric(selector) || is.logical(selector) || is.numeric(selector)) {
                selected_indices = selector
            }
            
            selected_indices
        },
        select_regions = function(selector) {
            self$get_regions()[private$select_region_indices(selector)]
        },
        select_region_metadata = function(selector) {
            private$region_metadata[private$select_region_indices(selector),, drop=F]
        }
    )
)
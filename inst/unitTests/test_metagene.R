## Test functions present in the metagene.R file

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "metagene2" )
    library( "data.table" )
    library( "dplyr" )
}

### }}}

bam_files <- get_demo_bam_files()
named_bam_files <- bam_files
names(named_bam_files) <- letters[1:(length(named_bam_files))]
not_indexed_bam_file <- metagene2:::get_not_indexed_bam_file()
regions <- metagene2:::get_demo_regions()
design <- data.frame(Samples = c("align1_rep1.bam", "align1_rep2.bam",
                    "align2_rep1.bam", "align2_rep2.bam", "ctrl.bam"),
                    align1 = c(1,1,0,0,2), align2 = c(0,0,1,1,2))
design$Samples <- paste0(system.file("extdata", package = "metagene2"), "/",
                        design$Samples)
set.seed(1)
demo_mg <- metagene2$new(regions = get_demo_regions(),
                        bam_files = get_demo_bam_files())
                        
full_mg = demo_mg$clone(deep=TRUE)
full_mg$produce_metagene()

test_invalid_param_constructor_value <- function(param_name, param_value, error_value) {
    arg_list_new = list(bam_files=get_demo_bam_files(), 
                        regions=get_demo_regions())
    arg_list_new[[param_name]] = param_value
    
    obs <- tryCatch(do.call(metagene2$new, arg_list_new),
                    error = conditionMessage)
    checkIdentical(obs, error_value)
}

test_invalid_param_value <- function(param_name, param_value, step_function, error_value) {
    arg_list_new = list(bam_files=get_demo_bam_files(), 
                        regions=get_demo_regions())
    arg_list_new[[param_name]] = param_value
    
    obs <- tryCatch(do.call(metagene2$new, arg_list_new),
                    error = conditionMessage)
    checkIdentical(obs, error_value)
    
    mg <- demo_mg$clone(deep=TRUE)
    obs <- tryCatch(do.call(mg[[step_function]], arg_list_new[param_name]),
                    error = conditionMessage)
    checkIdentical(obs, error_value)
    
    mg <- demo_mg$clone(deep=TRUE)
    obs <- tryCatch(do.call(mg$produce_metagene, arg_list_new[param_name]),
                    error = conditionMessage)
    checkIdentical(obs, error_value)      
}
                        
#region <- regions[1]
#bam_file <- bam_files[1]
#demo_mg_min <- metagene2$new(regions = region, bam_files = bam_file)

# Load fake SAMs for value tests.
fake_bam_files = system.file(c("fake_align1.bam", "fake_align2.bam", "fake_align3.bam"),
                             package="metagene2")
# Define a design that will test various combinations of bam files.
fake_bam_design = data.frame(BAM=fake_bam_files,
                             fake_align1=c(1, 0, 0),
                             fake_align2=c(0, 1, 0),
                             fake_align3=c(0, 0, 1),
                             fake_align12=c(1, 1, 0),
                             fake_align13=c(1, 0, 1),
                             fake_align23=c(0, 1, 1),
                             fake_align123=c(1, 1, 1),
                             stringsAsFactors=FALSE)

# Define the region over which single-region tests will be run.
# fake_align1 has 1 read over the whole region.
# fake_align2 has 2 reads over the whole region.
# fake_align3 has 4 reads over the first half of the region, 8 over the second half.
fake_bam_region_1_chr="chr1"
fake_bam_region_1_pos_start=1000000
fake_bam_region_1_pos_end = 1004999
fake_bam_region_1_pos_range=fake_bam_region_1_pos_start:fake_bam_region_1_pos_end
fake_bam_region_1_test_region_unique = GRanges(paste0(fake_bam_region_1_chr, ":", fake_bam_region_1_pos_start, "-", fake_bam_region_1_pos_end))

# Load the expected genomic coverages generated along with the bam files.
fake_bam_expected_coverages = list()
fake_bam_expected_rpm = list()
for(i in fake_bam_design$BAM) {
    cov_file = gsub(".bam", ".sam", i)
    load(paste0(cov_file, ".coverage.RData"))
    fake_bam_expected_coverages[[i]] = cov_rle
    
    load(paste0(cov_file, ".coverage_rpm.RData"))
    fake_bam_expected_rpm[[i]] = rpm_rle
}

###################################################
## Test the metagene2$new() function (initialize)
###################################################

## Invalid verbose value
test.metagene_invalid_verbose <- function() {
    test_invalid_param_constructor_value("verbose", "ZOMBIES", "verbose must be a logicial value (TRUE or FALSE)")
}

## Invalid force_seqlevels value
test.metagene_invalid_force_seqlevels <- function() {
    test_invalid_param_constructor_value("force_seqlevels", "ZOMBIES", "force_seqlevels must be a logicial value (TRUE or FALSE)")
}

## Negative padding_size value
test.metagene_invalid_padding <- function() {
    test_invalid_param_constructor_value("padding_size", "ZOMBIES", "padding_size must be a non-negative integer")
    test_invalid_param_constructor_value("padding_size", -1, "padding_size must be a non-negative integer")
    test_invalid_param_constructor_value("padding_size", 1.2, "padding_size must be a non-negative integer")
}


## Negative core value
test.metagene_invalid_core <- function() {
    test_invalid_param_constructor_value("cores", "ZOMBIES", "cores must be a positive integer or a BiocParallelParam instance.")
    test_invalid_param_constructor_value("cores", -1, "cores must be a positive integer or a BiocParallelParam instance.")
    test_invalid_param_constructor_value("cores", 1.2, "cores must be a positive integer or a BiocParallelParam instance.")
    test_invalid_param_constructor_value("cores", 0, "cores must be a positive integer or a BiocParallelParam instance.")    
}


## Non-character vector bam_files value
test_invalid_bam_file <- function(value, error_value) {
    obs <- tryCatch(metagene2:::metagene2$new(regions = get_demo_regions(),
                                            bam_files = value),
                    error = conditionMessage)
    checkIdentical(obs, error_value)
}

test.metagene_invalid_bam_files <- function() {
    test_invalid_bam_file(c(2,4,3), "bam_files must be a vector of BAM filenames.")
    test_invalid_bam_file(list(a = "a.txt", b = "b.txt"), "bam_files must be a vector of BAM filenames.")
    test_invalid_bam_file(not_indexed_bam_file, "All BAM files must be indexed")
    test_invalid_bam_file(c(bam_files, not_indexed_bam_file), "All BAM files must be indexed")
}

test_invalid_region <- function(invalid_region, error, ...) {
    obs <- tryCatch(do.call(metagene2$new, c(list(...), 
                                             list(bam_files=get_demo_bam_files(),
                                                  region = invalid_region))),
                    error = conditionMessage)
    checkIdentical(obs, error)
}

test_valid_region <- function(valid_region, ...) {
    obs <- do.call(metagene2$new, c(list(...), 
                                    list(bam_files=get_demo_bam_files(),
                                         region = valid_region)))
    test_valid_metagene(obs)
}

test_valid_metagene <- function(mg) {
    checkIdentical(class(mg), c("metagene", "R6"))
}

# regions is the wrong class
test.metagene_initialize_invalid_region_value <- function() {
    region <- array(data = NA, dim = c(2,2,2))
    test_invalid_region(region, paste0("regions must be either a vector of BED filenames, a ",
                                       "GRanges object or a GrangesList object"))
}

test.metagene_regions_seqlevels <- function() {
    # GRanges have seqlevels (indicating all possible seqnames)
    # and seqnames themselves. Extra seqlevels are no concerns.
    region_with_extra_seq_level = get_demo_regions()[[1]]
    GenomeInfoDb::seqlevels(region_with_extra_seq_level) <- c(GenomeInfoDb::seqlevels(region_with_extra_seq_level),
                                                              "extra_seqlevels")
    
    test_invalid_region(region_with_extra_seq_level, "Some seqlevels of regions are absent in bam_file")
    test_valid_region(region_with_extra_seq_level, force_seqlevels = TRUE)

    # Extra seqnames means we must drop certain regions. We only do so
    # if force_seqlevels=TRUE, otherwise we throw an error.
    region_with_extra_seq = region_with_extra_seq_level
    seqnames(region_with_extra_seq)[1] = "extra_seqlevels"

    test_invalid_region(region_with_extra_seq, "Some seqlevels of regions are absent in bam_file")
    test_valid_region(region_with_extra_seq, force_seqlevels = TRUE)

    # Sometimes there can be no seqnames left after removing
    # those with unknown levels. This happens often in chromosome names
    # are mismatched ("chr1" vs "1")
    region_no_common_seq = region_with_extra_seq_level
    seqnames(region_no_common_seq) = "extra_seqlevels"
    
    test_invalid_region(region_with_extra_seq, "No seqlevels matching between regions and bam file", force_seqlevels=TRUE)
}


# Valid regions narrowPeak
test.metagene_initialize_valid_narrowpeak <- function() {
    region <- metagene2:::get_narrowpeak_file()
    mg <- metagene2$new(regions = region, bam_files = bam_files[1])
    obs <- mg$get_regions()$list1
    extraCols <- c(signalValue = "numeric", pValue = "numeric",
                   qValue = "numeric", peak = "integer")
    exp <- rtracklayer::import(region, format = "BED", extraCols = extraCols)
    checkIdentical(obs, exp)
}

# Valid regions broadPeak
test.metagene_initialize_valid_broadpeak <- function() {
    region <- metagene2:::get_broadpeak_file()
    mg <- metagene2$new(regions = region, bam_files = bam_files[1])
    obs <- mg$get_regions()$list1
    extraCols <- c(signalValue = "numeric", pValue = "numeric",
                   qValue = "numeric")
    exp <- rtracklayer::import(region, format = "BED", extraCols = extraCols)
    checkIdentical(obs, exp)
}

# Valid named bam files
test.metagene_initialize_valid_bam_files <- function() {
    mg <- metagene2$new(regions = get_demo_regions(), bam_files = get_demo_bam_files())
    
    # Make sure all bam_files were kept.
    obs <- mg$get_params()[["bam_files"]]
    exp <- get_demo_bam_files()
    checkIdentical(obs, exp)
    
    # Make sure we have coverage for all bam files.
    obs <- names(mg$get_raw_coverages())
    exp <- get_demo_bam_files()
    checkIdentical(obs, exp)
    
    # NAmed bam files should have no impact (Historically, they changed the coverage names.)
    named_bam_files = get_demo_bam_files()
    names(named_bam_files) = letters[seq_along(named_bam_files)]
    mg <- metagene2$new(regions = get_demo_regions(), bam_files = named_bam_files)
    
    # Make sure we have coverage for all bam files under the correct names.
    obs <- names(mg$get_raw_coverages())
    exp <- unname(named_bam_files)
    checkIdentical(obs, exp)    
}

##################################################
# Test the metagene2$get_params() function
##################################################

## Valid usage
test.metagene_get_params_valid_usage <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    params <- mg$get_params()
    checkIdentical(unname(params[["bam_files"]]), get_demo_bam_files())
    checkIdentical(params[["padding_size"]], 0)
    checkIdentical(params[["verbose"]], FALSE)
    checkIdentical(params[["force_seqlevels"]], FALSE)
}

##################################################
# Test the metagene2$get_design() function
##################################################

## Valid usage
test.metagene_get_design_valid_usage <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    mg$group_coverages(design=get_demo_design())
    design <- mg$get_design()
    checkIdentical(mg$get_design()[,-1], get_demo_design()[,-1])
    checkIdentical(mg$get_design()[,1], mg$get_params()[["bam_files"]])
}

##################################################
# Test the metagene2$get_regions() function
##################################################

## Valid usage default
test.metagene_get_regions_valid_usage_default <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    regions <- mg$get_regions()
    exp <- get_demo_regions()
    exp <- tools::file_path_sans_ext(basename(exp))
    checkIdentical(length(regions), length(unlist(get_demo_regions)))
    region_1_in_regions = length(get_demo_regions()[[1]])
    region_1_in_metagene = sum(regions$region_name==names(get_demo_regions())[1])
    checkIdentical(region_1_in_regions, region_1_in_metagene)
}

##################################################
# Test the metagene2$get_matrice() function
##################################################

# Tests to write:
#  Test single region
#  Test bin_size = 1
#  Test bin_count < width(regions)

test.metagene_group_coverages_valid_usage_default = function(){
    mg <- demo_mg$clone(deep=TRUE)
    grouped_coverages = mg$group_coverages(design=get_demo_design())

    checkIdentical(length(grouped_coverages), ncol(get_demo_design()) - 1)
}

# Make sure bam_files order does not change results.
# See issue #7: https://github.com/ArnaudDroitLab/metagene2/issues/7
test.metagene_group_coverages_bam_files_order = function(){
    mg  = metagene2$new(regions=get_demo_regions(), bam_files=get_demo_bam_files())
    mg2 = metagene2$new(regions=get_demo_regions(), bam_files=rev(get_demo_bam_files()))
    checkIdentical(mg$bin_coverages()[["align1_rep1"]][1:10,1:10],
                   mg2$bin_coverages()[["align1_rep1"]][1:10,1:10])
}


test.metagene_bin_coverages_valid_usage_default_design = function(){
    mg <- demo_mg$clone(deep=TRUE)
    mg$group_coverages()
    binned_coverages = mg$bin_coverages()

    design = mg$get_design()
    checkIdentical(length(binned_coverages), ncol(design) - 1)
    checkIdentical(names(binned_coverages), colnames(design)[-1])
    for(i in 1:length(binned_coverages)) {
        checkIdentical(nrow(binned_coverages[[i]]) == length(mg$get_regions()))
        checkIdentical(ncol(binned_coverages[[i]]) == mg$get_params()[["bin_count"]])
    }
}

test.metagene_bin_coverages_valid_usage_custom_design = function(){
    mg <- demo_mg$clone(deep=TRUE)
    mg$group_coverages(design=get_demo_design())
    binned_coverages = mg$bin_coverages()

    checkIdentical(length(binned_coverages), ncol(get_demo_design()) - 1)
    checkIdentical(names(binned_coverages), colnames(get_demo_design())[-1])
    for(i in 1:length(binned_coverages)) {
        checkIdentical(nrow(binned_coverages[[i]]) == length(mg$get_regions()))
        checkIdentical(ncol(binned_coverages[[i]]) == mg$get_params()[["bin_count"]])
    }
}

test.metagene_split_coverages_valid_usage_default = function(){
    mg <- demo_mg$clone(deep=TRUE)
    mg$group_coverages(design=get_demo_design())
    mg$split_coverages()

    checkIdentical(length(binned_coverages), ncol(get_demo_design()) - 1)
    checkIdentical(names(binned_coverages), colnames(get_demo_design())[-1])
    for(i in 1:length(binned_coverages)) {
        checkIdentical(nrow(binned_coverages[[i]]) == length(mg$get_regions()))
        checkIdentical(ncol(binned_coverages[[i]]) == mg$get_params()[["bin_count"]])
    }
}

test.metagene_raw_coverage_values_unique_region = function(){
    # Create the metagene object.
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, regions=test_region_unique)
    
    # Test genome-wide raw and normalized coverages for
    # strand_mode=FALSE
    raw_coverages = mg$get_raw_coverages()
    norm_coverages = mg$get_normalized_coverages()
    for(i in fake_bam_design$BAM) {
        obs = raw_coverages[[i]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
        exp = fake_bam_expected_coverages[[i]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
        checkIdentical(obs, exp)
                       
        obs = norm_coverages[[i]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
        exp = fake_bam_expected_rpm[[i]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
        checkTrue(all(abs(obs - exp) < 10e-8))
    }
}

test.metagene_group_coverage_values_unique_region = function(){
    # Create the metagene object.
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, regions=test_region_unique)
     
    # Test group coverages.
    group_coverages = mg$group_coverages(design=fake_bam_design)
    
    # Grouped coverages for unit designs should be the same as the plain ones.
    for(i in c("fake_align1", "fake_align2", "fake_align3")) {
        obs = group_coverages[["*"]][[i]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
        exp = fake_bam_expected_coverages[[paste0(i, ".bam")]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
        checkIdentical(obs, exp)
    }

    # Test combined group coverages.
    obs = group_coverages[["*"]][["fake_align12"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    exp = fake_bam_expected_coverages[["fake_align1.bam"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range] +
          fake_bam_expected_coverages[["fake_align2.bam"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    checkIdentical(obs, exp)
    
    obs = group_coverages[["*"]][["fake_align23"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    exp = fake_bam_expected_coverages[["fake_align2.bam"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range] +
          fake_bam_expected_coverages[["fake_align3.bam"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    checkIdentical(obs, exp)
    
    obs = group_coverages[["*"]][["fake_align13"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    exp = fake_bam_expected_coverages[["fake_align1.bam"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range] +
          fake_bam_expected_coverages[["fake_align3.bam"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    checkIdentical(obs, exp)

    obs = group_coverages[["*"]][["fake_align123"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    exp = fake_bam_expected_coverages[["fake_align1.bam"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range] +
          fake_bam_expected_coverages[["fake_align2.bam"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range] +
          fake_bam_expected_coverages[["fake_align3.bam"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    checkIdentical(obs, exp)
}

test.metagene_bin_coverage_values_unique_region = function(){
    # Create the metagene object.
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, regions=test_region_unique, design=fake_bam_design)
    
    # Test binned coverages for bins 
    # where all values are the same.
    bin_coverages = mg$bin_coverages()
    checkTrue(all(bin_coverages[["fake_align1"]][1,]==1))
    checkTrue(all(bin_coverages[["fake_align2"]][1,]==2))
    checkTrue(all(bin_coverages[["fake_align3"]][1,1:50]==4))
    checkTrue(all(bin_coverages[["fake_align3"]][1,51:100]==8))
    checkTrue(all(bin_coverages[["fake_align123"]][1,1:50]==7))
    checkTrue(all(bin_coverages[["fake_align123"]][1,51:100]==11))
    
    # Test binned coverages where some values are means.
    bin_coverages = mg$bin_coverages(bin_count=5)
    checkTrue(all(bin_coverages[["fake_align1"]][1,]==1))
    checkTrue(all(bin_coverages[["fake_align2"]][1,]==2))
    checkTrue(bin_coverages[["fake_align3"]][1,1]==4)
    checkTrue(bin_coverages[["fake_align3"]][1,2]==4)
    checkTrue(bin_coverages[["fake_align3"]][1,3]==6)
    checkTrue(bin_coverages[["fake_align3"]][1,4]==8)
    checkTrue(bin_coverages[["fake_align3"]][1,5]==8)
}

test.metagene_calculate_ci_values_unique_region = function(){
    # Create the metagene object.
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, regions=test_region_unique, design=fake_bam_design)
 
    full_df = mg$add_metadata()

    # Only one region, all confidence intervals should be NA.
    checkTrue(all(is.na(full_df$qinf) & is.na(full_df$qsup)))
    
    # There should be five bins.
    checkTrue(max(full_df$bin)==5 && length(unique(full_df$bin))==5)
    
    # Bin values should match
    checkTrue((full_df %>% filter(bin==3 & design=="fake_align3") %>% pull(value)) == 6)
    
    # All groups should be represented.
    checkTrue(all(full_df$design %in% colnames(fake_bam_design)[-1]))
    
    test_plot = mg$produce_metagene()
    checkTrue("ggproto" %in% class(test_plot$scales))
}


##################################################
# Test the metagene2$get_data_frame() function
##################################################

## Valid usage default
test.metagene_get_data_frame_valid_usage_default <- function() {
 df <- full_mg$get_data_frame()
 regions <- get_demo_regions()
 bam_files <- get_demo_bam_files()
 checkTrue(is.data.frame(df))
 checkTrue(ncol(df) == 8)
 checkTrue(nrow(df) == length(regions) * length(bam_files) * mg$get_params()$bin_count)
}

### Valid usage subset
test.metagene_get_data_frame_valid_usage_subset <- function() {
 single_region = names(get_demo_regions())[1]
 three_designs = colnames(full_mg$get_design())[2:4]

 df <- full_mg$get_data_frame(region_names = single_region, design_names =three_designs )

 checkTrue(is.data.frame(df))
 checkTrue(ncol(df) == 8)
 checkTrue(nrow(df) == length(single_region) * length(three_designs) * full_mg$get_params()$bin_count)
}

## Valid usage get_data_frame return by copy of data_frame
test.metagene_get_data_frame_check_copy_of_data_frame <- function() {
 mg <- full_mg$clone(deep=TRUE)
 
 df1 <- mg$get_data_frame()
 #modification of table by reference
 df1$c <- 1:1000

 #Is table copied and unchanged ? 
 df2 <- mg$get_data_frame()
 checkIdentical(ncol(df1) == ncol(df2), FALSE)
}

## Valid usage no data_frame produced
test.metagene_get_data_frame_valid_usage_no_data_frame <- function() {
 mg <- demo_mg$clone(deep=TRUE)
 df <- mg$get_data_frame()
 checkTrue(is.null(df))
 df_subset <- mg$get_data_frame(get_demo_regions()[1],
                                get_demo_bam_files()[1:2])
 checkTrue(is.null(df_subset))
}

##################################################
# Test the metagene2$get_plot() function
##################################################

## Valid case no graph
test.metagene_get_plot_valid_case_no_graph <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    plot <- mg$get_plot()
    checkTrue(is.null(plot))
}

## Valid case graph
test.metagene_get_plot_valid_case_graph <- function() {
    plot <- full_mg$get_plot()
    checkTrue(all(class(plot) == c("gg", "ggplot")))
}

##################################################
# Test the metagene2$get_raw_coverages() function
##################################################

exp_raw <- GenomicAlignments::readGAlignments(bam_files[1])
exp_raw <- GenomicAlignments::coverage(exp_raw)

## Default filenames
test.metagene_get_raw_coverages_default_filenames <- function() {
    obs <- full_mg$get_raw_coverages()[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) all(obs[[i]]==exp_raw[[i]]),
                        logical(1))))
}

## One filename
test.metagene_get_raw_coverages_one_filename <- function() {
    obs <- full_mg$get_raw_coverages(filenames = bam_files[1])[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) all(obs[[i]]==exp_raw[[i]]),
                        logical(1))))
}

## All filenames
test.metagene_get_raw_coverages_all_filename <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- mg$get_raw_coverages(filenames = bam_files)[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) all(obs[[i]]==exp_raw[[i]]),
                        logical(1))))
}

## Invalid filenames class
test.metagene_get_raw_coverages_invalid_filenames_class <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- tryCatch(mg$get_raw_coverages(filenames = 1),
                    error = conditionMessage)
    exp <- "is.character(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid empty filename
test.metagene_get_raw_coverages_invalid_empty_filename <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- tryCatch(mg$get_raw_coverages(filenames = ""),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid filename alone
test.metagene_get_raw_coverages_invalid_filename_alone <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- tryCatch(mg$get_raw_coverages(filenames = "asdf"),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid filename among valid
test.metagene_get_raw_coverages_invalid_filename_among_valid <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- tryCatch(mg$get_raw_coverages(filenames = c("asdf", bam_files)),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

##################################################
# Test the metagene2$get_normalized_coverages() function
##################################################

count <- Rsamtools::countBam(bam_files[1])$records
weight <- 1 / (count / 1000000)
exp_norm <- exp_raw * weight

## Default filenames
test.metagene_get_normalized_coverages_default_filenames <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- mg$get_normalized_coverages()[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) all(obs[[i]]==exp_raw[[i]]),
                        logical(1))))
}

## One filename
test.metagene_get_normalized_coverages_one_filename <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- mg$get_normalized_coverages(filenames = bam_files[1])[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) all(obs[[i]]==exp_raw[[i]]),
                        logical(1))))
}

## All filenames
test.metagene_get_normalized_coverages_all_filename <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- mg$get_normalized_coverages(filenames = bam_files)[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) all(obs[[i]]==exp_raw[[i]]),
                        logical(1))))
}

## Invalid filenames class
test.metagene_get_normalized_coverages_invalid_filenames_class <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- tryCatch(mg$get_normalized_coverages(filenames = 1),
                    error = conditionMessage)
    exp <- "is.character(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid empty filename
test.metagene_get_normalized_coverages_invalid_empty_filename <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- tryCatch(mg$get_normalized_coverages(filenames = ""),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid filename alone
test.metagene_get_normalized_coverages_invalid_filename_alone <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- tryCatch(mg$get_normalized_coverages(filenames = "asdf"),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid filename among valid
test.metagene_get_normalized_coverages_invalid_filename_among_valid <- 
    function() {
    mg <- demo_mg$clone(deep=TRUE)
    filenames <- c("asdf", bam_files)
    obs <- tryCatch(mg$get_normalized_coverages(filenames = filenames),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

##################################################
# Test the metagene2$add_design() function
##################################################

## Valid design data frame
test.metagene_add_design_valid_design_data_frame <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    mg$add_design(get_demo_design())
    checkIdentical(mg$get_design()[,-1], get_demo_design()[,-1])
    checkIdentical(mg$get_design()[,1], names(mg$get_params()[["bam_files"]]))
}

## Valid design NULL
test.metagene_add_design_valid_design_null <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    mg$add_design(design = NULL)
    checkIdentical(colnames(mg$get_design())[-1],
                names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid design NA, NA first
test.metagene_add_design_valid_design_na_na_first <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    mg$add_design(design = NA)
    checkIdentical(colnames(mg$get_design())[-1],
                names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid design NA, NULL first
test.metagene_add_design_valid_design_na_null_first <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    mg$add_design(design = NULL)
    mg$add_design(design = NA)
    checkIdentical(colnames(mg$get_design())[-1],
                names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid design NA, design first
test.metagene_add_design_valid_design_na_design_first <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    mg$add_design(design = get_demo_design())
    checkIdentical(mg$get_design()[,-1], get_demo_design()[,-1])
    checkIdentical(mg$get_design()[,1], names(mg$get_params()[["bam_files"]]))
}

## Valid design, factor sample names
test.metagene_add_design_valid_design_factor_sample_names <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    design <- get_demo_design()
    design[,1] <- factor(design[,1])
    mg$add_design(design)
    checkTrue(is.factor(design[,1]))
    checkTrue(is.character(mg$get_design()[,1]))
    checkIdentical(design[,-1], mg$get_design()[,-1])
    checkIdentical(names(mg$get_params()[["bam_files"]]), mg$get_design()[,1])
}

## Valid check_bam_files TRUE NA design
test.metagene_add_design_valid_check_bam_files_true_na_design <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    mg$add_design(design = NA, check_bam_files = TRUE)
    checkIdentical(colnames(mg$get_design())[-1],
                names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid check_bam_files TRUE NULL design
test.metagene_add_design_valid_check_bam_files_true_null_design <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    mg$add_design(design = NA, check_bam_files = TRUE)
    checkIdentical(colnames(mg$get_design())[-1],
                names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid check_bam_files TRUE design design
test.metagene_add_design_valid_check_bam_files_true_design_design <- function()
{
    mg <- demo_mg$clone(deep=TRUE)
    mg$add_design(design = get_demo_design(), check_bam_files = TRUE)
    checkIdentical(mg$get_design()[,-1], get_demo_design()[,-1])
    checkIdentical(mg$get_design()[,1], names(mg$get_params()[["bam_files"]]))
}

## Invalid design class
test.metagene_add_design_invalid_design_class <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- tryCatch(mg$add_design(design = 1), error = conditionMessage)
    exp <- "design must be a data.frame object, NULL or NA"
    checkIdentical(obs, exp)
}

## Invalid design column
test.metagene_add_design_invalid_design_column <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    design <- get_demo_design()
    design <- design[, 1, drop = FALSE]
    obs <- tryCatch(mg$add_design(design = design), error = conditionMessage)
    exp <- "design must have at least 2 columns"
    checkIdentical(obs, exp)
}

## Invalid design column one class
test.metagene_add_design_invalid_design_column_one_class <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    design <- get_demo_design()
    design[,1] <- seq_along(design[,1])
    obs <- tryCatch(mg$add_design(design = design), error = conditionMessage)
    exp <- "The first column of design must be BAM filenames"
    checkIdentical(obs, exp)
}

## Invalid design column two plus class
test.metagene_add_design_invalid_design_columns_two_plus_class <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    design <- get_demo_design()
    design[,2] <- letters[seq_along(design[,2])]
    obs <- tryCatch(mg$add_design(design = design), error = conditionMessage)
    exp <- "All design column, except the first one, must be in numeric format"
    checkIdentical(obs, exp)
}

## Invalid bam file check_bam_files TRUE
test.metagene_add_design_invalid_bam_file_check_bam_files_true <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    design <- get_demo_design()
    design[1,1] <- "not_a_valid_bam_file"
    obs <- tryCatch(mg$add_design(design = design, check_bam_files = TRUE),
                    error = conditionMessage)
    exp <- "Design contains samples absent from the list of bam files provided on initialization."
    checkIdentical(obs, exp)
}

## Invalid bam file check_bam_files FALSE
test.metagene_add_design_invalid_bam_file_check_bam_files_false <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    design <- get_demo_design()
    design[1,1] <- "not_a_valid_bam_file"
    obs <- tryCatch(mg$add_design(design = design, check_bam_files = TRUE),
                    error = conditionMessage)
    exp <- "Design contains samples absent from the list of bam files provided on initialization."
    checkIdentical(obs, exp)
}

##################################################
# Test the metagene2$produce_table() function
##################################################

test.metagene_produce_table_valid_without_design <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    checkIdentical("bin_count" %in% mg$get_params(), FALSE)
    mg$produce_table()
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(is.data.frame(mg$get_table()), TRUE)
    #length of table : number of region * number of design * number of bin * number of range by region (demo = 50,000 lines)
    tablength <- length(mg$get_regions())*length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    checkIdentical(dim(mg$get_table())[1] == tablength, TRUE)
    checkIdentical(dim(mg$get_table())[2] == length(c('region', 'design', 'bin', 'value', 'strand')), TRUE)
    tab <- mg$get_table()
    #check for presence of levels of factors (region, design, strand)
    checkIdentical(names(mg$get_regions()), unique(tab$region))
    checkIdentical(tools::file_path_sans_ext(basename(mg$get_params()$bam_files)), unique(tab$design))
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(tab$strand[which(tab$region == region_names)]))
    }
    #check for number of line by factor region
    reglength <- length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(tab$region[which(tab$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in tools::file_path_sans_ext(basename(mg$get_params()$bam_files))){
            #print(design_names)
            checkIdentical(length(tab$design[which(tab$region == region_names & tab$design == design_names)]) == designlength , TRUE)
        }
    }
    print(TRUE)
}

test.metagene_produce_table_valid_with_design <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    checkIdentical("bin_count" %in% mg$get_params(), FALSE)
    demo_design <- get_demo_design()
    mg$produce_table(design = demo_design)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(is.data.frame(mg$get_table()), TRUE)
    #length of table : number of region * number of design * number of bin * number of range by region (demo = 50,000 lines)
    tablength <- length(mg$get_regions())*(dim(demo_design)[2]-1)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    checkIdentical(dim(mg$get_table())[1] == tablength, TRUE)
    checkIdentical(dim(mg$get_table())[2] == length(c('region', 'design', 'bin', 'value', 'strand')), TRUE)
    tab <- mg$get_table()
    #check for presence of levels of factors (region, design, strand)
    checkIdentical(names(mg$get_regions()), unique(tab$region))
    checkIdentical(names(demo_design)[-1], unique(tab$design))
    
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(tab$strand[which(tab$region == region_names)]))
    }
    #check for number of line by factor region
    reglength <- (dim(demo_design)[2]-1)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(tab$region[which(tab$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in names(demo_design)[-1]){
            #print(design_names)
            checkIdentical(length(tab$design[which(tab$region == region_names & tab$design == design_names)]) == designlength , TRUE)
        }
    }
    print(TRUE)
}

test.metagene_produce_table_valid_without_design_bin_count_50 <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    checkIdentical("bin_count" %in% mg$get_params(), FALSE)
    mg$produce_table(bin_count = 50)
    checkIdentical(mg$get_params()[["bin_count"]], 50)
    checkIdentical(is.data.frame(mg$get_table()), TRUE)
    #length of table : number of region * number of design * number of bin * number of range by region (demo = 50,000 lines)
    tablength <- length(mg$get_regions())*length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    checkIdentical(dim(mg$get_table())[1] == tablength, TRUE)
    checkIdentical(dim(mg$get_table())[2] == length(c('region', 'design', 'bin', 'value', 'strand')), TRUE)
    tab <- mg$get_table()
    #check for presence of levels of factors (region, design, strand)
    checkIdentical(names(mg$get_regions()), unique(tab$region))
    checkIdentical(tools::file_path_sans_ext(basename(mg$get_params()$bam_files)), unique(tab$design))
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(tab$strand[which(tab$region == region_names)]))
    }
    #check for number of line by factor region
    reglength <- length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(tab$region[which(tab$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in tools::file_path_sans_ext(basename(mg$get_params()$bam_files))){
            #print(design_names)
            checkIdentical(length(tab$design[which(tab$region == region_names & tab$design == design_names)]) == designlength , TRUE)
        }
    }
    print(TRUE)
}

test.metagene_produce_table_valid_with_design_bin_count_50 <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    checkIdentical("bin_count" %in% mg$get_params(), FALSE)
    demo_design <- get_demo_design()
    mg$produce_table(design = demo_design, bin_count = 50)
    checkIdentical(mg$get_params()[["bin_count"]], 50)
    checkIdentical(is.data.frame(mg$get_table()), TRUE)
    #length of table : number of region * number of design * number of bin * number of range by region (demo = 50,000 lines)
    tablength <- length(mg$get_regions())*(dim(demo_design)[2]-1)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    checkIdentical(dim(mg$get_table())[1] == tablength, TRUE)
    checkIdentical(dim(mg$get_table())[2] == length(c('region', 'design', 'bin', 'value', 'strand')), TRUE)
    tab <- mg$get_table()
    #check for presence of levels of factors (region, design, strand)
    checkIdentical(names(mg$get_regions()), unique(tab$region))
    checkIdentical(names(demo_design)[-1], unique(tab$design))
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(tab$strand[which(tab$region == region_names)]))
    }
    #check for number of line by factor region
    reglength <- (dim(demo_design)[2]-1)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(tab$region[which(tab$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in names(demo_design)[-1]){
            #print(design_names)
            checkIdentical(length(tab$design[which(tab$region == region_names & tab$design == design_names)]) == designlength , TRUE)
        }
    }
    print(TRUE)
}
#
# Not valid design object
test.metagene_produce_table_invalid_design <- function() {
 mg <- demo_mg$clone(deep=TRUE)
 obs <- tryCatch(mg$produce_table(design = c(1,2)),
                error = conditionMessage)
 exp <- "design must be a data.frame object, NULL or NA"
 checkIdentical(obs, exp)
}
#
# Design data.frame with not enough columns
test.metagene_produce_table_invalid_design_data_frame <- function() {
 mg <- demo_mg$clone(deep=TRUE)
 design <- data.frame(a = c("ZOMBIE_ONE", "ZOMBIE_TWO"))
 obs <- tryCatch(mg$produce_table(design = design),
                error = conditionMessage)
 exp <- "design must have at least 2 columns"
 checkIdentical(obs, exp)
}

# Design data.frame with invalid first column
test.metagene_produce_table_invalid_design_first_column <- function() {
 mg <- demo_mg$clone(deep=TRUE)
 design <- data.frame(a = c(1,3), zombies = c("ZOMBIE_ONE", "ZOMBIE_TWO"))
 obs <- tryCatch(mg$produce_table(design = design),
                error = conditionMessage)
 exp <- "The first column of design must be BAM filenames"
 checkIdentical(obs, exp)
}

# Design data.frame with invalid second column
test.metagene_produce_table_invalid_design_second_column <- function() {
 mg <- demo_mg$clone(deep=TRUE)
 designTemp<-data.frame(a = named_bam_files,
                        zombies = rep("ZOMBIE_ONE", length(named_bam_files)))
 obs <- tryCatch(mg$produce_table(design = designTemp),
                error = conditionMessage)
 exp <- paste0("All design column, except the first one, must be in ",
                "numeric format")
 checkIdentical(obs, exp)
}

# Design data.frame with invalid second column
test.metagene_produce_table_invalid_design_not_defined_file <- function() {
 mg <- demo_mg$clone(deep=TRUE)
 designNew<-data.frame(a = c(bam_files, "I am not a file"),
                        b = rep(1, length(bam_files) + 1))
 obs <- tryCatch(mg$produce_table(design = designNew),
                error = conditionMessage)
 exp <- "Design contains samples absent from the list of bam files provided on initialization."
 checkIdentical(obs, exp)
}

# Design using zero file (0 in all rows of the design object)
test.metagene_produce_table_design_using_no_file <- function() {
 mg <- demo_mg$clone(deep=TRUE)
 designNew<-data.frame(a = bam_files,
                        b = rep(0, length(bam_files)))
 obs <- tryCatch(mg$produce_table(design = designNew),
                error = conditionMessage)
 exp <- "At least one BAM file must be used in the design."
 checkIdentical(obs, exp)
}

# Invalid bin_count class
test.metagene_invalid_bin_count <- function() {
  test_invalid_param_value("bin_count", "a", "bin_coverages", "bin_count must be a positive integer")
  test_invalid_param_value("bin_count", -1, "bin_coverages", "bin_count must be a positive integer")
  test_invalid_param_value("bin_count", 1.2, "bin_coverages", "bin_count must be a positive integer")  
}



# Valid noise_removal NCIS
test.metagene_produce_table_valid_noise_removal_ncis <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    design <- get_demo_design()[,1:2]
    design[,2][2] <- 0
    mg$produce_table(noise_removal = "NCIS", design = design)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(mg$get_params()[["noise_removal"]], "NCIS")
    tab <- mg$get_table()
    tablength <- length(mg$get_regions())*(dim(design)[2]-1)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    checkIdentical(dim(tab)[1] == tablength, TRUE)
    checkIdentical(dim(tab)[2] == length(c('region', 'design', 'bin', 'value', 'strand')), TRUE)
    
    tab <- mg$get_table()
    #check for presence of levels of factors (region, design, strand)
    checkIdentical(names(mg$get_regions()), unique(tab$region))
    checkIdentical(names(design)[-1], unique(tab$design))
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(tab$strand[which(tab$region == region_names)]))
    }
    #check for number of line by factor region
    reglength <- (dim(design)[2]-1)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(tab$region[which(tab$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in names(design)[-1]){
            #print(design_names)
            checkIdentical(length(tab$design[which(tab$region == region_names & tab$design == design_names)]) == designlength , TRUE)
        }
    }
    print(TRUE)
}

# Invalid normalization class
test.metagene_invalid_normalization <- function() {
 test_invalid_param_value("normalization", 1234, "group_coverages", 'normalization must be NULL, "RPM" or "NCIS".')
 test_invalid_param_value("normalization", "CSI", "group_coverages", 'normalization must be NULL, "RPM" or "NCIS".')
}

## Valid normalization RPM
test.metagene_produce_table_valid_normalization_rpm <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    mg$produce_table(normalization = "RPM")
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(mg$get_params()[["normalization"]], "RPM")
    tab <- mg$get_table()
    tablength <- length(mg$get_regions())*length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    checkIdentical(dim(tab)[1] == tablength, TRUE)
    checkIdentical(dim(tab)[2] == length(c('region', 'design', 'bin', 'value', 'strand')), TRUE)
    
    tab <- mg$get_table()
    #check for presence of levels of factors (region, design, strand)
    checkIdentical(names(mg$get_regions()), unique(tab$region))
    checkIdentical(tools::file_path_sans_ext(basename(mg$get_params()$bam_files)), unique(tab$design))
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(tab$strand[which(tab$region == region_names)]))
    }
    #check for number of line by factor region
    reglength <- length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(tab$region[which(tab$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in tools::file_path_sans_ext(basename(mg$get_params()$bam_files))){
            #print(design_names)
            checkIdentical(length(tab$design[which(tab$region == region_names & tab$design == design_names)]) == designlength , TRUE)
        }
    }
    print(TRUE)
}
#
## Invalid flip_regions class
test.metagene_produce_table_invalid_flip_regions_class <- function() {
 mg <- demo_mg$clone(deep=TRUE)
 obs <- tryCatch(mg$produce_table(flip_regions = 1234),
                error = conditionMessage)
 exp <- "flip_regions must be a logical."
 checkIdentical(obs, exp)
}
#
## Valid flip_regions true
test.metagene_produce_table_valid_flip_regions_true <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$produce_table(flip_regions = TRUE)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
    
    #modifier strand char for regions
    
    
    #test expected == observed
    tab <- mg$get_table()
    # print(tab)
    expect <- c()
    for (region_names in names(mg$get_regions())){
        if(as.vector(strand(mg$get_regions()[[region_names]]))[1] == "-") {
            expect <- c(expect,rep(100:1,length(names(mg$get_design())[-1])*length(as.vector(strand(mg$get_regions()[[region_names]])))))
        } else {
            expect <- c(expect,rep(1:100,length(names(mg$get_design())[-1])*length(as.vector(strand(mg$get_regions()[[region_names]])))))
        }
    }
    # print(class(tab$bin))
    # print(class(expect))
    checkIdentical(as.numeric(tab$bin), as.numeric(expect))
}
#
## Valid flip_regions false
test.metagene_produce_table_valid_flip_regions_false <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$produce_table(flip_regions = FALSE)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    
    #modifier strand char dans regions
    
    
    #test expected == observed
    tab <- mg$get_table()
    # print(tab)
    expect <- c()
    for (region_names in names(mg$get_regions())){
        # if(as.vector(strand(mg$get_regions()[[region_names]]))[1] == "-") {
            # expect <- c(expect,rep(100:1,length(names(mg$get_design())[-1])*length(as.vector(strand(mg$get_regions()[[region_names]])))))
        # } else {
            expect <- c(expect,rep(1:100,length(names(mg$get_design())[-1])*length(as.vector(strand(mg$get_regions()[[region_names]])))))
        # }
    }
    # print(class(tab$bin))
    # print(class(expect))
    checkIdentical(as.numeric(tab$bin), as.numeric(expect))
}

##################################################
# Test the metagene2$produce_data_frame() function
##################################################

## Valid default usage
test.metagene_produce_data_frame_default_arguments <- function(){
    mg <- demo_mg$clone(deep=TRUE)
    mg$produce_table()
    mg$produce_data_frame()
    df <- mg$get_data_frame()
    #check the nrow & ncol
    checkIdentical(ncol(df), length(c('region', 'design', 'bin', 'value', 'strand', 'qinf', 'qsup', 'group')))
    expectedNbRow <- length(mg$get_regions())*length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)
    checkIdentical(nrow(df) == expectedNbRow, TRUE)
    #check colnames
    checkIdentical(colnames(df), c('region', 'design', 'bin', 'value', 'strand', 'qinf', 'qsup','group'))
    #check region, design repartition in data_frame
    checkIdentical(names(mg$get_regions()), unique(df$region))
    checkIdentical(names(mg$get_design()[-1]), unique(df$design))
    #check for number of line by factor region
    reglength <- length(names(mg$get_design())[-1])*(mg$get_params()$bin_count)
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(df$region[which(df$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in names(mg$get_design()[-1])){
            #print(design_names)
            checkIdentical(length(df$design[which(df$region == region_names & df$design == design_names)]) == designlength , TRUE)
        }
    }
    #check for bin repartition
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(sum(df$bin), sum(1:100)*length(names(mg$get_design())[-1])*length(mg$get_regions()))
    }
    #strand matches
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(df$strand[which(df$region == region_names)]))
    }
    print(TRUE)
}


### Invalid alpha value
test.metagene_invalid_alpha <- function(){
    test_invalid_param_value("alpha", 'test', "calculate_ci", "alpha >= 0 & alpha <= 1 is not TRUE")
    test_invalid_param_value("alpha", -0.8, "calculate_ci", "alpha >= 0 & alpha <= 1 is not TRUE")    
}

### Invalid sample_count class
test.metagene_invalid_sample_count <- function(){
    test_invalid_param_value("sample_count", 'test', "calculate_ci", "is.numeric(sample_count) is not TRUE")
    test_invalid_param_value("sample_count", -10, "calculate_ci", "sample_count >= 0 is not TRUE")    
}

test.metagene_replace_metadata <- function() {
    regions_gr <- unlist(get_demo_regions())
    demo_metadata = data.frame(BedName=names(regions_gr),
                               EvenStart=ifelse((start(regions_gr) %% 2) == 0, "Even", "Odd"),
                               Strand=strand(regions_gr))

    mg <- metagene2$new(regions = get_demo_regions(),
                       region_metadata=demo_metadata,
                       bam_files = get_demo_bam_files(),
                       assay='chipseq')
    
    mg$produce_metagene(split_by=c("EvenStart", "Strand"))
    test_meta = mg$get_regions_metadata()
    
    # Add column and replace. This should work fine.
    test_meta$Foo = c("Foo", "Bar")
    mg$replace_region_metadata(test_meta)

    # Remove region_name column. It should be added back.
    old_region_name = mg$get_regions_metadata()$region_name
    test_meta = test_meta[setdiff(colnames(test_meta), "region_name")]
    mg$replace_region_metadata(test_meta)
    
    checkTrue(all(mg$get_regions_metadata()$region_name==old_region_name))
    
    # Remove one of the split_by columns. Caches should be invalidated.
    test_meta = test_meta[setdiff(colnames(test_meta), "EvenStart")]
    mg$replace_region_metadata(test_meta)    
    checkIdentical(mg$get_params()[["split_by"]], "region_name")
    
    # Try using a data-frame without enough rows. We should get an error message.
    obs <- tryCatch(mg$replace_region_metadata(test_meta[1:4,]),
                    error = conditionMessage)
    exp <- "region_metadata must have one row per region."
    checkIdentical(obs, exp)    
    
}
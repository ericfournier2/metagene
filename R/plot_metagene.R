#' Produce a metagene plot
#' 
#' @param df a \code{data.frame} obtained with the \code{get_data_frame}
#' function. Must have the following columns: "region", "design", "bin", 
#' "value", "qinf" and "qsup".
#'
#' @return A `ggplot` object.
#'
#' @examples
#' region <- get_demo_regions()[1]
#' bam_file <- get_demo_bam_files()[1]
#' mg <- metagene$new(regions = region, bam_files = bam_file)
#' mg$produce_data_frame()
#' df <- mg$get_data_frame()
#' p <- plot_metagene(df)
plot_metagene <- function(df, facet_by=NULL, group_by=NULL) {
    df$design <- as.factor(df$design)

    expected_cols <- c("bin", "value", "qinf", "qsup", "group")
    assert_subset<-df[,which(colnames(df) %in% expected_cols)]
    expected_class <- c("integer", rep("numeric", 3), "factor")
    names(expected_class) = expected_cols
    stopifnot(all(expected_cols %in% colnames(assert_subset)))
    actual_classes = vapply(assert_subset, class, character(1))
    actual_classes = actual_classes[expected_cols]
    stopifnot(all(actual_classes == expected_class))

    if(is.null(group_by)) {
        group_by="group"
    }
    
    p <- ggplot(df, aes(x=bin, y=value, ymin=qinf, ymax=qsup)) +
        geom_ribbon(aes_string(fill=group_by), alpha=0.3) +
        geom_line(aes_string(color=group_by), size=1) +
        theme(panel.grid.major = element_line()) +
        theme(panel.grid.minor = element_line()) +
        theme(panel.background = element_blank()) +
        theme(panel.background = element_rect()) +
        theme_bw(base_size = 20)
        
    if(!is.null(facet_by)) {
        p <- p + facet_grid(facet_by)
    }
    
    return(p)
}
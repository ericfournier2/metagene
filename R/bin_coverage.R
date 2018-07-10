tile_indices <- function(run_length, n, width) {
    if(!is.null(n)) {
        return(tile(IRanges(start=1, end=run_length), n=n))
    } else {
        return(tile(IRanges(start=1, end=run_length), width=width))
    }
}

tile_indices_start <- function(run_length, n, width) {
    unlist(start(tile_indices(run_length, n, width)))
}

tile_indices_end <- function(run_length, n, width) {
    unlist(end(tile_indices(run_length, n, width)))
}

flip_regions_list <- function(list, cond) {
    mapply(list, cond, FUN=function(x,y) {
      if(y) {
          return(rev(x))
      } else {
          return(x)
      }    
  })    
}

reduce_sum_group <- function(x, group_length) {
    stopifnot((length(x) %% group_length) == 0)
    total = rep(0, group_length)
    for(i in seq(1, (length(x) / group_length))) {
        index_base = ((i - 1) * group_length) + 1
        total <- total + x[index_base:((index_base+group_length) - 1)]
    }
    return(total)
}

#bin_uniform_contiguous_regions <- function(coverage, regions, bin_count) {
#  region_views = Views(coverage, start(regions), end(regions))
#  
#  region_strands = as.character(strand(regions))
#  region_rles = flip_regions_list(region_views, region_strands=="-")
#  bin_rles = lapply(region_rles, function(x) {
#    unlist(lapply(split(x, rep(1:bin_count, each=length(x) / bin_count)), sum))
#  })
#  res_matrix = matrix(unlist(bin_rles), ncol=100, byrow=TRUE)
#  
#  region_sum = Reduce("+", region_rles)
#  tile_starts = tile_indices_start(length(region_sum), n=bin_count)
#  tile_ends = tile_indices_end(length(region_sum), n=bin_count)
#  bin_views = Views(region_sum, start=tile_starts, end=tile_ends)
#  
#  return(list(Matrix=res_matrix / (width(regions[1]) / bin_count),
#              Vector=viewMeans(bin_views) / length(regions)))
#}

# Harder case: Continuous ranges of heteregeneous lengths
bin_heterogeneous_contiguous_regions <- function(coverage, regions, bin_count) {
  region_strands = as.character(strand(regions))
  tile_regions_list = tile(regions, n=bin_count)
  tile_regions = unlist(GRangesList(tile_regions_list))
  tile_views = Views(coverage, start(tile_regions), end(tile_regions))
  tile_means = viewMeans(tile_views)
  m <- matrix(tile_means, ncol = bin_count, byrow = TRUE)
  mr <- m[,bin_count:1]
  i <- region_strands == "-"
  m[i,] <- mr[i,]
#  bin_means <- colSums(m) / length(regions)
#  return(list(Matrix=m, Vector=bin_means))
  m
}

# Hardest case: Discontinuous ranges of heteregeneous lengths
bin_heterogeneous_discontiguous_regions <- function(coverage, regions_list, bin_count) {
    region_rles = lapply(regions_list, function(x) {
        unlist(Views(coverage, start(x), end(x)))
    })
    region_strands = lapply(regions_list, function(x) {
        return(as.character(strand(x)[1]))
    })
    region_rles = mapply(region_rles, region_strands, FUN=function(x,y) {
        if(y=="-") {
            return(rev(x))
        } else {
            return(x)
        }    
    })
    region_rles_views = lapply(region_rles, function(x) {
        Views(x, tile_indices_start(length(x), n=bin_count), tile_indices_end(length(x), n=bin_count))
    })
    region_bin_means = lapply(region_rles_views, viewMeans)
    bin_means = reduce_sum_group(unlist(region_bin_means), bin_count)  / length(regions_list)
    return(list(Matrix=matrix(unlist(region_bin_means), ncol=bin_count, byrow=TRUE),
                Vector=bin_means))
    
}

get_subtable = function(coverages, gr, bcount) {
    grl <- split(gr, GenomeInfoDb::seqnames(gr))
    i <- vapply(grl, length, numeric(1)) > 0
    do.call("c", lapply(grl[i], get_view_means,
                        bcount = bcount, cov = coverages))
}
get_view_means = function(gr, bcount, cov) {
    chr <- unique(as.character(GenomeInfoDb::seqnames(gr)))
    gr <- unlist(tile(unlist(gr), n=bcount))
    stopifnot(length(chr) == 1)
    views <- Views(cov[[chr]], start(gr), end(gr))
    viewMeans(views)
}

bin_contiguous_regions <- function(coverage, regions, bin_count) {
  m <-  matrix(get_subtable(coverage, regions, bin_count), ncol=bin_count, byrow=TRUE)
    
  mr <- m[,bin_count:1]
  i <-as.logical(strand(unlist(regions_gr))=="-")
  m[i,] <- mr[i,]
  
  m
}

bin_coverages = function(coverages, regions, bin_count) {
    results = list()
    for(cov_name in names(coverages)) {
        if(is(regions, "GRangesList")) {
            results[[cov_name]] = bin_heterogeneous_discontiguous_regions(coverages[[cov_name]], regions, bin_count)
        } else {
            results[[cov_name]] = bin_contiguous_regions(coverages[[cov_name]], regions, bin_count)
        }        
    }
    return(results)
}

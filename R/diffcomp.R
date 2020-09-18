#' Calculate differential compartments
#'
#' This function uses biological replicates (bioreps) to calculate a
#' false-discovery rate (fdr) for determining differential compartments between
#' two samples.
#'
#' Pair-wise euclidean distances between input vectors are calculated between
#' sample-matched biological replicates (bioreps) and biorep-matched samples.
#' All distance vectors are then combined and ranked by decreasing value. The
#' false-discovery rate (fdr) is calculated cumulatively along the ranked
#' distance vector as a percentage of bioreps. Significance is determined by
#' finding the greatest rank that passes the fdr cutoff threshold.
#'
#' @param data data.frame with columns of compartment scores and rows for each
#'   Hi-C bin.
#'
#' @param colData data.frame describing the columns of \code{data}. It must
#'   contain columns describing samples and bioreps. \code{rownames(colData)}
#'   should equal \code{colnames(data)}.
#'
#' @param featureData (optional) data.frame describing the rows of \code{data}.
#'   Each column of \code{featureData} should describe different features of the
#'   corresponding rows in \code{data}. \code{nrow(featureData)} should equal
#'   \code{nrow(data)}.
#'
#' @param bioreps string naming a column in \code{colData} that describes which
#'   rows belong to which biological replicate.
#'
#' @param samples string naming a column in \code{colData} that describes which
#'   rows belong to which sample.
#'
#' @param cutoff numeric threshold for false-discovery rate.
#'
#' @return Returns two data.frames. \code{distdf} contains \code{featureData},
#'   \code{data}, and the pair-wise distances between bioreps and samples along
#'   with a column denoting significance (T/F). \code{ranked} contains a ranked
#'   listing of samples, bioreps, fdr, and rank.
#'
#' @examples
#' ## Load example compartment data
#' data("compartments")
#'
#' ## Define eigenvector data
#' eigen <- compartments[,grep("EI", colnames(compartments))]
#'
#' ## Define colData describing samples and bioreps
#' cData <- data.frame(
#'   time = c("0000", "0000", "1440", "1440"),
#'   br = c(1, 2, 1, 2),
#'   row.names = colnames(eigen)
#' )
#'
#' ## Define featureData describing rows of data
#' fData <- compartments[,1:3]
#'
#' ## Calculate differential compartments
#' res <- diffcomp(data = eigen,
#'                 colData = cData,
#'                 featureData = fData,
#'                 bioreps = "br",
#'                 samples = "time",
#'                 cutoff = 0.01)
#'
#' ## Look at results
#' lapply(res, head)
#'
#' ## Filter for significant compartments
#' sigComp <- res$distdf[res$distdf$significant == T,]
#' head(sigComp)
#' nrow(sigComp)
#'
#'
#'
#' @export
diffcomp <- function(data, colData, featureData, bioreps, samples, cutoff = 0.05) {

  ## Extract factor levels from colData
  biorep_levels <- levels(as.factor(colData[[bioreps]]))
  sample_levels <- levels(as.factor(colData[[samples]]))

  ## Calculate biorep distances between matched samples
  bio <- lapply(1:length(sample_levels), function(i) {
    apply(data[rownames(colData)[colData[samples] == sample_levels[i]]], 1, dist)
  })

  ## Calculate sample distances between matched bioreps
  sam <- lapply(1:length(biorep_levels), function(i) {
    apply(data[rownames(colData)[colData[bioreps] == biorep_levels[i]]], 1, dist)
  })

  ## Add names to each set of distances
  names(bio) <- paste0("bio_", paste0(samples, sample_levels), "_")
  names(sam) <- paste0("sam_", paste0(bioreps, biorep_levels), "_")

  ## Assemble into a data frame
  df <- data.frame(
    do.call(cbind, bio),
    do.call(cbind, sam)
  )

  ## Rank distances from highest to lowest; rename columns; add rank
  ranked <- as.data.frame(
    do.call(rbind, strsplit(names(sort(unlist(df), decreasing = T)), "_"))
  )

  ## Calculate false-discovery rate for each cumulative rank
  fdr <- cumsum(ranked[,1] == "bio") / 1:nrow(ranked)

  ## Add fdr, rank and name columns
  ranked[,4] <- fdr
  ranked[,5] <- 1:nrow(ranked)
  colnames(ranked) <- c("comparison", "constant", "distdf_row", "fdr", "rank")

  ## Find the maximum rank that meets the cutoff
  sigRank <- max(which(ranked$fdr <= cutoff))

  ## Annotate significant regions
  df$significant <- F
  df$significant[unique(ranked$distdf_row[ranked$rank <= sigRank])] <- T

  ## Format output data (distdf) ##
  ## Fix column names
  colnames(df) <- gsub("_$", "", colnames(df))

  ## Add featureData (if supplied)
  if(!missing(featureData)) distdf <- featureData

  ## Add original data and distance data
  distdf <- data.frame(distdf, data, df)

  ## Return results
  return(list(
    distdf = distdf,
    ranked = ranked
  ))
}




## Bug check that rownames of colData match colnames(data)

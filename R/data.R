#' Compartment Data
#'
#' A dataset containing compartment calls for
#' chromosome 1 at 100Kb resolution using
#' Eigenvector from Juicer tools. Compartments
#' were called in 2 timepoints (i.e. 0000 and 1440)
#' as well as 2 biological replicates
#' (i.e. BR2 and BR3) of a Hi-C timecourse for
#' macrophages activated with LPS + INF-g.
#'
#' @format An object of class data.frame with 2493 rows and 7 columns.
#' \describe{
#'   \item{seqnames}{chromosome name}
#'   \item{start}{starting position of each 100Kb bin}
#'   \item{end}{ending position of each 100Kb bin}
#'   \item{EI_0000_BR2}{eigenvector value for time = 0000 and biorep 2}
#'   \item{EI_0000_BR3}{eigenvector value for time = 0000 and biorep 3}
#'   \item{EI_1440_BR2}{eigenvector value for time = 1440 and biorep 2}
#'   \item{EI_1440_BR3}{eigenvector value for time = 1440 and biorep 3}
#' }
#'
#' @source Eigenvector from Juicer tools:
#'  \url{https://github.com/aidenlab/juicer/wiki/Eigenvector}
#'
"compartments"

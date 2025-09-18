#' Read Highland Tribes
#'
#' A static network of political alliances and enmities among the 16 Gahuku-Gama
#' sub-tribes of Eastern Central Highlands of New Guinea, documented by \insertCite{read1954cultures;textual}{ergm.sign}.
#'
#' @name tribes
#' @docType data
#' @format
#'  An undirected \code{static.sign} object with no loops.
#' @references Taken from UCINET IV, which cites the following: \insertRef{hage1983structural}{ergm.sign}  \insertRef{read1954cultures}{ergm.sign}
#' @source
#' \url{http://vlado.fmf.uni-lj.si/pub/networks/data/UciNet/UciData.htm#gama},
#' with corrections from \insertCite{read1954cultures;textual}{ergm.sign}.
#' @examples
#'
#' \donttest{
#' data(tribes)
#' }
#'
NULL

#' Conflict Events in Syrian Civil War
#'
#' A dynamic network of combat events in the Syrian civil war between 2017 and
#' 2019. The dataset is taken from \insertCite{fritz2023all;textual}{ergm.sign}.
#' The raw data for this example Armed Conflict Location & Event Data Project \insertCite{raleigh2010introducing;textual}{ergm.sign}.
#'
#' @name rebels
#' @docType data
#' @format An undirected \code{dynamic.sign} object with no loops and eight timepoints.
#' @references \insertRef{fritz2023all}{ergm.sign} \insertRef{raleigh2010introducing}{ergm.sign}
#'
#' @examples
#'
#' \donttest{
#' data(rebels)
#' }
NULL


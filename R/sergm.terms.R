#' Terms Used in Signed Exponential Family Random Graph Models
#'
#' @name sergm.terms
#' @description How to specify netowrk statistics in the [`ergm.sign`][ergm.sign] package.
#' @section Specifying models:
#' The lhs of the formula needs to specify an object of the class \code{static.sign} or \code{dynamic.sign}.
#' Similarly to the definition of terms in \CRANpkg{ergm}, the rhs contains a description of the networks
#' terms in an additive manner. As for the naming convention of the currently implemented
#' network terms, the suffix \eqn{_pos} indicates that the network term only regard positive edges,
#' while the suffix \eqn{_neg} does the same for negative edges. If there is none of these two suffices
#' in the name of a term, it relates to positive and negative edges at the same time. In addition to the usual \link{ergm} terms, some terms have been added specifically for analyzing signed networks in relation to structural balance theory.
#'
#' @section Edges:
#' \enumerate{
#' \item \code{edges_pos}:
#'
#' This adds a term counting all positive edges, i.e., where the adjacency matrix network equals 1,
#'  in the specified network.
#' \item \code{edges_neg}:
#'
#' This adds a term counting all negative edges, i.e., where the adjacency matrix network equals -1,
#'  in the specified network.
#' \item \code{edges}:
#'
#' This adds a term counting any type of edges, i.e., where the adjacency matrix network is unequal to 0,
#'  in the specified network.
#' }
#'
#' @section Isolates:
#' \enumerate{
#' \item \code{isolates_pos}:
#'
#' This adds a term counting all actors in the network with no positive edges.
#' \item \code{isolates_neg}:
#'
#' This adds a term counting all actors in the network with no negative edges.
#' \item \code{isolates}:
#'
#' This adds a term counting all actors in the network with no edges, be they positive or negative.
#' }
#'
#' @section Degree:
#' \enumerate{
#' \item \code{degree_pos(d = c(i:j))}:
#'
#' This adds a separate term counting all actors in the networks that have positive degree of i, i+1, ..., j-1, and j.
#'
#' \item \code{degree_neg(d = c(i:j))}:
#'
#' This adds a separate term counting all actors in the networks that have negative degree of i, i+1, ..., j-1, and j.
#'
#' \item \code{degree(d = c(i:j))}:
#'
#' This adds a separate term counting all actors in the networks that have degree of i, i+1, ..., j-1, and j.
#'
#' \item \code{gwdegree_pos(decay = alpha, fixed = FALSE, attrname = NULL, cutoff = 30)}:
#'
#' This adds a term to the model of the geometrically weighted positive degrees.
#' The statistic is equal the sum of actors with a specific number of positive degree,
#' and the number of actors with degree k is weighted by exp(alpha)*(1-(1-exp(-alpha))^k).
#
#' \item \code{gwdegree_neg(decay = alpha, fixed = FALSE, attrname = NULL, cutoff = 30)}:
#'
#' This adds a term to the model of the geometrically weighted negative degrees.
#' The statistic is equal the sum of actors with a specific number of negative degree,
#' and the number of actors with degree k is weighted by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwdegree(decay = alpha, fixed = FALSE, attrname = NULL, cutoff = 30)}:
#'
#' This adds a term to the model of the geometrically weighted degrees.
#' The statistic is equal the sum of actors with a specific number of degree,
#' and the number of actors with degree k is weighted by exp(alpha)*(1-(1-exp(-alpha))^k).
#' The degree of actor n is defined as the number of positive and negative ties
#' actor n has in the network.
#' }
#'
#' @section Edgewise-shared partners:
#' \enumerate{
#'
#'\item \code{esf_pos(d, type="OTP")}:
#' This adds a term for the count of positive edgewise-shared friends to the model.
#' The count is the number of friends that each have a positive tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#' Relating to the structural balance theory this term translates to
#' clustering according to the 'friends-of-friends-are-friends' mechanism.
#'
#' \item \code{esf_neg(d, type="OTP")}:
#' This adds a term for the count of negative edgewise-shared friends to the model.
#' The count is the number of enemies that each have a positive tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#' Relating to the structural balance theory this term translates to
#' clustering according to the 'friends-of-enemies-are-friends' mechanism.
#'
#' \item \code{ese_pos(d, type="OTP")}:
#' This adds a term for the count of positive edgewise-shared enemies to the model.
#' The count is the number of friends that each have a negative tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#' Relating to the structural balance theory this term translates to
#' clustering according to the 'enemies-of-friends-are-enemies' mechanism.
#'
#' \item \code{ese_neg(d, type="OTP")}:
#' This adds a term for the count of negative edgewise-shared enemies to the model.
#' The count is the number of enemies that each have a negative tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#' Relating to the structural balance theory this term translates to
#' clustering according to the 'enemies-of-enemies-are-enemies' mechanism.
#'
#' #' \item \code{esm_pos(d, type="OTP")}:
#' This adds a term for the count of positive edgewise-shared mixed partners to the model.
#' The count is the number of friends that have a positive and a negative tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#' Relating to the structural balance theory this term translates to
#' clustering according to the 'friends-of-friends-are-enemies' mechanism.
#'
#' \item \code{esm_neg(d, type="OTP")}:
#' This adds a term for the count of negative edgewise-shared mixed partners to the model.
#' The count is the number of enemies that have a positive and a negative tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#' Relating to the structural balance theory this term translates to
#' clustering according to the 'friends-of-enemies-are-enemies' mechanism.
#'
#'  \item \code{gwesf_pos(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted positive edgewise-shared friends.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwesf_neg(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted negative edgewise-shared friends.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwese_pos(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted positive edgewise-shared enemies.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwese_neg(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted negative edgewise-shared enemies.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwesm_pos(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted positive edgewise-shared mixed partners.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwesm_neg(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted negative edgewise-shared mixed partners.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#' }
#'
#' @section Dyadwise-shared partners:
#' \enumerate{
#' \item \code{dsf(d, type="OTP")}:
#' This adds a term for the count of dyadwise-shared friends to the model.
#' The count is the number of pairs of actors (connected or unconnected) that each have a positive tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#'
#' \item \code{dse(d, type="OTP")}:
#' This adds a term for the count of dyadwise-shared enemies to the model.
#' The count is the number of pairs of actors (connected or unconnected) that each have a negative tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#'
#' \item \code{dsm(d, type="OTP")}:
#' This adds a term for the count of dyadwise-shared mixed partners to the model.
#' The count is the number of pairs of actors (connected or unconnected) that have a positive and a negative tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#'
#' \item \code{gwdsf(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted dyadwise-shared friends
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwdse(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted dyadwise-shared enemies.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwdsm(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted dyadwise-shared mixed partners.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#' }
#'
#' @section Non-edgewise-shared partners:
#' \enumerate{
#' \item \code{nesf_pos(d, type="OTP")}:
#' This adds a term for the count of positive non-edgewise-shared friends to the model.
#' The count is the number of pairs of actors that are not friends but each have a positive tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#'
#' \item \code{nesf_neg(d, type="OTP")}:
#' This adds a term for the count of negative non-edgewise-shared friends to the model.
#' The count is the number of pairs of actors that are not enemies but each have a positive tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#'
#' \item \code{nese_pos(d, type="OTP")}:
#' This adds a term for the count of positive non-edgewise-shared enemies to the model.
#' The count is the number of pairs of actors that are not friends but each have a negative tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#'
#' \item \code{nese_neg(d, type="OTP")}:
#' This adds a term for the count of negative non-edgewise-shared enemies to the model.
#' The count is the number of pairs of actors that are not enemies but each have a negative tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#'
#' \item \code{nesm_pos(d, type="OTP")}:
#' This adds a term for the count of positive non-edgewise-shared mixed partners to the model.
#' The count is the number of pairs of actors that are not friends but are connected through a positive and negative tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#'
#' \item \code{nesm_neg(d, type="OTP")}:
#' This adds a term for the count of negative non-edgewise-shared mixed partners to the model.
#' The count is the number of pairs of actors that are not enemies but are connected through a positive and negative tie to a common third actor, this value can lie between 0 and
#' n -2 (where n is the number of actors in the network).
#'
#'  \item \code{gwnesf_pos(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted positive non-edgewise-shared friends.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwnesf_neg(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted negative non-edgewise-shared friends.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwnese_pos(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted positive non-edgewise-shared enemies.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwnese_neg(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted negative non-edgewise-shared enemies.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwnesm_pos(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted positive non-edgewise-shared mixed partners.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#'
#' \item \code{gwnesm_neg(decay = alpha, fixed = FALSE, cutoff = 30, type = "OTP")}:
#' This adds a term to the model of the geometrically weighted negative non-edgewise-shared mixed partners.
#' The weight is given by exp(alpha)*(1-(1-exp(-alpha))^k).
#' }
#'
NULL

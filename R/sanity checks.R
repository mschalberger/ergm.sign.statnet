#' Consistency Checks for Package
#'
#' Performs consistency checks for package
#'
#'@name sanity
#'
#' @examples
#' set.seed(1234)
#' mat <- matrix(sample(c(-1,0,1),100,replace=TRUE),ncol=10)
#' mat[lower.tri(mat)] <- 0
#' diag(mat) <- 0
#' static_net <- signnet(mat, directed = F, loops = F, matrix.type = "adjacency")
#'
#' # Check that uncombine_network gets the same solution for Layer as manually
#' set.seed(123)
#' model <- sergm(static_net ~ edges_pos + edges_neg + gwesf_pos(decay=0.5, fix = T))
#'
#' # uncombine_network
#' uncomb <- uncombine_network(simulate(model, seed = 123), split.vattr = ".LayerName")
#' neg <- as.sociomatrix(uncomb[[1]])*-1
#' pos <- as.sociomatrix(uncomb[[2]])
#' comb <- neg+pos
#'
#' # manually
#' net <- as.matrix.network(simulate(model, seed = 123))
#' n <- ncol(net)
#' comb2 <- net[1:(n/2),1:(n/2)] + net[(n/2+1):n,(n/2+1):n]*-1
#'
#' stopifnot(identical(comb, comb2))
#'
#' # gwespL and gwesp on Layer is the same
#' class(static_net) <- "network"
#' MultiLayer <- Layer(static_net, c(`+` = "pos",`-`= "neg"))
#' set.seed(123)
#' ergm(MultiLayer ~ gwespL(0.5, fixed = T, Ls.path = c(~`+`,~`+`), L.base = ~ `+`))
#' summary_formula(MultiLayer ~ gwespL(0.5, fixed = T, Ls.path = c(~`+`,~`+`), L.base = ~ `+`))
#' set.seed(123)
#' ergm(MultiLayer ~ L(~ gwesp(0.5, fixed = T) , Ls = ~`+`))
#' summary_formula(MultiLayer ~ L(~ gwesp(0.5, fixed = T) , Ls = ~`+`))
#'
#' # impact of L.in_order
#' set.seed(123)
#' ergm(MultiLayer ~ gwespL(0.5, fixed = T, Ls.path = c(~`+`,~`-`), L.base = ~ `+`, L.in_order = T))
#' set.seed(123)
#' ergm(MultiLayer ~ gwespL(0.5, fixed = T, Ls.path = c(~`+`,~`-`), L.base = ~ `+`, L.in_order = F))
#'

NULL


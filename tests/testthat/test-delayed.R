library(testthat)
library(ergm)

context("Delayed ERGM terms with multiple previous networks via summary(dynamic.sign)")

make_multi_networks <- function() {
  # --- t1 ---
  adj1 <- matrix(0, 5, 5)
  adj1[1,3] <- 1
  adj1[2,3] <- 1
  adj1[3,4] <- -1
  adj1[3,5] <- -1
  adj1[lower.tri(adj1)] <- t(adj1)[lower.tri(adj1)]
  net1 <- network.sign(adj1, matrix.type = "adjacency", directed = F)
  net1%v%"party" <- c("A", "A", "B", "B", "B")

  # --- t2 ---
  adj2 <- matrix(0, 5, 5)
  adj2[1,2] <- -1
  adj2[2,3] <- -1
  adj2[3,4] <- 1
  adj2[4,5] <- 1
  adj2[lower.tri(adj1)] <- t(adj2)[lower.tri(adj2)]
  net2 <- network.sign(adj2, matrix.type = "adjacency", directed = F)
  net2%v%"party" <- c("A", "B", "A", "B", "B")

  # --- t3 ---
  adj3 <- matrix(0, 5, 5)
  adj3[1,3] <- -1
  #adj3[1,2] <- -1
  adj3[4,5] <- 1
  adj3[3,5] <- 1
  adj3[lower.tri(adj3)] <- t(adj3)[lower.tri(adj3)]
  net3 <- network.sign(adj3, matrix.type = "adjacency", directed = F)
  net3%v%"party" <- c("A", "A", "B", "B")

  nets <- networks.sign(list(net1,net2,net3), dynamic = T)
  class(nets) <-class(nets)[-1]

  return(nets)
}

make_multi_networks_directed <- function() {
  # --- t1 ---
  adj1 <- matrix(0, 5, 5)
  adj1[1,3] <- 1
  adj1[2,3] <- 1
  adj1[3,4] <- -1
  net1 <- network.sign(adj1, matrix.type = "adjacency", directed = TRUE)

  # --- t2 ---
  adj2 <- matrix(0, 5, 5)
  adj2[3,1] <- 1
  adj2[3,2] <- 1
  adj2[4,3] <- -1
  net2 <- network.sign(adj2, matrix.type = "adjacency", directed = TRUE)

  # --- t3 ---
  adj3 <- matrix(0, 5, 5)
  adj3[1,3] <- 1
  adj3[4,5] <- 1
  adj3[5,3] <- 1
  adj3[3,4] <- -1
  net3 <- network.sign(adj3, matrix.type = "adjacency", directed = TRUE)
  net3%v%"party" <- c("A", "A", "B", "B", "B")

  nets <- networks.sign(list(net1, net2, net3), dynamic = TRUE)
  class(nets) <-class(nets)[-1]
  return(nets)
}

# --- Tests -------------------------------------------------------------

test_that("delayed ese and and gwese work as expected and give the same value for delay = 0 and no specified d", {
  net <- make_multi_networks()
  s_pos <- summary(net ~ Cross(~delese(base="+")))
  s_neg <- summary(net ~ Cross(~delese(base="-")))
  gw_s_pos <- summary(net ~ Cross(~gwdelese(decay=0, base="+")))
  gw_s_neg <- summary(net ~ Cross(~gwdelese(decay=0, base="-")))
  expect_equal(unname(gw_s_pos), 1)
  expect_equal(unname(s_pos), 1)
  expect_equal(unname(gw_s_neg), 1)
  expect_equal(unname(s_neg), 1)
})

test_that("delayed esf and and gwes work as expected and give the same value for delay = 0 and no specified d", {
  net <- make_multi_networks()
  gw_s_pos <- summary(net ~ Cross(~delesf(base="+")))
  gw_s_neg <- summary(net ~ Cross(~delesf(base="-")))
  s_pos <- summary(net ~ Cross(~gwdelesf(decay=0, base="+")))
  s_neg <- summary(net ~ Cross(~gwdelesf(decay=0, base="-")))
  expect_equal(unname(gw_s_pos), 1)
  expect_equal(unname(s_pos), 1)
  expect_equal(unname(gw_s_neg), 1)
  expect_equal(unname(s_neg), 1)
})


test_that("delnodematch works for all layers", {
  net <- make_multi_networks()
  all <- summary(net ~ Cross(~delnodematch("party")))
  pos <- summary(net ~Cross(~Pos(~delnodematch("party"))))
  neg <- summary(net ~Cross(~Neg(~delnodematch("party"))))
  expect_equal(unname(all), unname(pos + neg))
})


test_that("delrecip works for all layers (directed)", {
  net <- make_multi_networks_directed()
  s_pos <- summary(net ~ Cross(~Pos(~delrecip)))
  s_neg <- summary(net ~ Cross(~Neg(~delrecip)))
  s_all <- summary(net ~ Cross(~delrecip))
  expect_equal(unname(s_all), unname(s_pos + s_neg))
})

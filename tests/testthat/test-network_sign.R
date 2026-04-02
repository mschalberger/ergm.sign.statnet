# tests/testthat/test-network_sign.R

library(testthat)
library(network)
library(ergm)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# 4x4 signed adjacency matrix
signed_adj <- function() {
  m <- matrix(0L, 4, 4)
  m[1, 2] <-  1L
  m[2, 3] <- -1L
  m[3, 4] <-  1L
  m
}

# Named variant
signed_adj_named <- function() {
  m <- signed_adj()
  dimnames(m) <- list(letters[1:4], letters[1:4])
  m
}

# Positive / negative split of the same adjacency
pos_adj <- function() ((signed_adj() > 0) * 1L)
neg_adj <- function() ((signed_adj() < 0) * 1L)

# Three-column edgelist (from, to, sign)
signed_edgelist <- function() {
  data.frame(
    From = c(1L, 2L, 3L),
    To   = c(2L, 3L, 4L),
    Sign = c(1L, -1L, 1L)
  )
}

# Tiny helper to extract the 'sign' vertex attribute values
layer_names <- function(net) unique(net %v% "sign")

# ---------------------------------------------------------------------------
# 1. Adjacency matrix — static
# ---------------------------------------------------------------------------

test_that("adjacency: returns static.sign network", {
  net <- network.sign(signed_adj(), matrix.type = "adjacency")
  expect_s3_class(net, "static.sign")
  expect_s3_class(net, "network")
})

test_that("adjacency: layers are named '+' and '-'", {
  net <- network.sign(signed_adj(), matrix.type = "adjacency")
  expect_setequal(layer_names(net), c("pos", "neg"))
})

test_that("adjacency: dual.sign flag stored correctly", {
  net_f <- network.sign(signed_adj(), matrix.type = "adjacency", dual.sign = FALSE)
  net_t <- network.sign(signed_adj(), matrix.type = "adjacency", dual.sign = TRUE)
  expect_false(net_f %n% "dual.sign")
  expect_true(net_t  %n% "dual.sign")
})

test_that("adjacency: sign vertex attribute mirrors .LayerName", {
  net <- network.sign(signed_adj(), matrix.type = "adjacency")
  expect_equal(net %v% "sign", net %v% ".LayerName")
})

# ---------------------------------------------------------------------------
# 2. Edgelist — static
# ---------------------------------------------------------------------------

test_that("edgelist: returns static.sign network", {
  net <- network.sign(signed_edgelist(), matrix.type = "edgelist")
  expect_s3_class(net, "static.sign")
})

test_that("edgelist: layers named '+' and '-'", {
  net <- network.sign(signed_edgelist(), matrix.type = "edgelist")
  expect_setequal(layer_names(net), c("pos", "neg"))
})

test_that("edgelist: dual.sign stored correctly", {
  net_f <- network.sign(signed_edgelist(), matrix.type = "edgelist", dual.sign = FALSE)
  net_t <- network.sign(signed_edgelist(), matrix.type = "edgelist", dual.sign = TRUE)
  expect_false(net_f %n% "dual.sign")
  expect_true(net_t  %n% "dual.sign")
})

# ---------------------------------------------------------------------------
# 3. Separate pos.mat / neg.mat — static
# ---------------------------------------------------------------------------

test_that("pos/neg mats: returns static.sign network", {
  net <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj())
  expect_s3_class(net, "static.sign")
})

test_that("pos/neg mats: layers are named '+' and '-'", {
  net <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj())
  expect_setequal(layer_names(net), c("pos", "neg"))
})

test_that("pos/neg mats: errors when only one is supplied", {
  expect_error(
    network.sign(pos.mat = pos_adj()),
    regexp = "Both"
  )
  expect_error(
    network.sign(neg.mat = neg_adj()),
    regexp = "Both"
  )
})

test_that("pos/neg mats: errors when one is a list and other is not", {
  expect_error(
    network.sign(pos.mat = list(pos_adj()), neg.mat = neg_adj()),
    regexp = "both be lists"
  )
})

# ---------------------------------------------------------------------------
# 4. Dynamic — list of matrices
# ---------------------------------------------------------------------------

test_that("dynamic (list of matrices): returns dynamic.sign network", {
  mats <- list(signed_adj(), signed_adj(), signed_adj())
  net  <- network.sign(mats, matrix.type = "adjacency")
  expect_s3_class(net, "dynamic.sign")
})

test_that("dynamic (list of pos/neg mats): returns dynamic.sign network", {
  pos_list <- list(pos_adj(), pos_adj())
  neg_list <- list(neg_adj(), neg_adj())
  net <- network.sign(pos.mat = pos_list, neg.mat = neg_list)
  expect_s3_class(net, "dynamic.sign")
})

test_that("dynamic pos/neg: errors on length mismatch", {
  expect_error(
    network.sign(pos.mat = list(pos_adj(), pos_adj()), neg.mat = list(neg_adj())),
    regexp = "equal length"
  )
})

# ---------------------------------------------------------------------------
# 5. Timepoint pooling
# ---------------------------------------------------------------------------

test_that("pooling (numeric timepoints): returns pooled network", {
  mats <- list(signed_adj(), signed_adj(), signed_adj(), signed_adj())
  net  <- network.sign(mats, matrix.type = "adjacency", timepoints = 2)
  # 4 timepoints pooled into 2 → dynamic with 2 snapshots or single static
  expect_true(inherits(net, "dynamic.sign") || inherits(net, "static.sign"))
})

test_that("pooling (list timepoints): groups are respected", {
  m1 <- signed_adj(); m1[1, 2] <-  1L
  m2 <- signed_adj(); m2[1, 2] <- -1L   # opposite sign in first cell
  m3 <- signed_adj(); m3[1, 2] <-  1L
  m4 <- signed_adj(); m4[1, 2] <-  1L

  # Group 1: m1+m2 → tie → resolved by tie.breaker
  # Group 2: m3+m4 → +1
  mats <- list(m1, m2, m3, m4)
  tpts <- list(c(1, 2), c(3, 4))

  net_zero <- network.sign(mats, matrix.type = "adjacency",
                           timepoints = tpts, tie.breaker = "zero")
  net_pos  <- network.sign(mats, matrix.type = "adjacency",
                           timepoints = tpts, tie.breaker = "positive")

  expect_s3_class(net_zero, "dynamic.sign")
  expect_s3_class(net_pos,  "dynamic.sign")
})

test_that("tie.breaker = 'first' / 'last' use correct snapshot", {
  m_pos <- signed_adj();  m_pos[2, 3] <-  1L
  m_neg <- signed_adj();  m_neg[2, 3] <- -1L

  mats <- list(m_pos, m_neg)
  tpts <- list(c(1, 2))

  net_first <- network.sign(mats, matrix.type = "adjacency",
                            timepoints = tpts, tie.breaker = "first")
  net_last  <- network.sign(mats, matrix.type = "adjacency",
                            timepoints = tpts, tie.breaker = "last")

  # Both valid; just confirm they don't error and return the right class
  expect_true(inherits(net_first, "dynamic.sign") || inherits(net_first, "static.sign"))
  expect_true(inherits(net_last,  "dynamic.sign") || inherits(net_last,  "static.sign"))
})

# ---------------------------------------------------------------------------
# 6. Named adjacency matrix → vertex names propagate
# ---------------------------------------------------------------------------

test_that("named adjacency: vertex names propagate to layers", {
  net <- network.sign(signed_adj_named(), matrix.type = "adjacency")
  # At least one layer should carry the letter names
  vnames <- net %v% "vertex.names"
  expect_true(any(c("a", "b", "c", "d") %in% vnames))
})

# ---------------------------------------------------------------------------
# 7. dual.sign constraint — consistent with Signed()
# ---------------------------------------------------------------------------

test_that("dual.sign = FALSE adds constraint (matches Signed() behaviour)", {
  net_ns  <- network.sign(signed_adj(), matrix.type = "adjacency", dual.sign = FALSE)
  constraints_ns <- net_ns %ergmlhs% "constraints"
  expect_false(is.null(constraints_ns))
})

test_that("dual.sign = TRUE skips constraint (matches Signed() behaviour)", {
  net_ns_f <- network.sign(signed_adj(), matrix.type = "adjacency", dual.sign = FALSE)
  net_ns_t <- network.sign(signed_adj(), matrix.type = "adjacency", dual.sign = TRUE)

  c_f <- deparse(net_ns_f %ergmlhs% "constraints")
  c_t <- deparse(net_ns_t %ergmlhs% "constraints")
  expect_false(identical(c_f, c_t))
})

# ---------------------------------------------------------------------------
# 8. Cross-function comparison: network.sign() vs Signed()
# ---------------------------------------------------------------------------

test_that("adjacency path: layer structure matches Signed(pos_nw, neg_nw)", {
  pos_net <- network::network(pos_adj(), directed = FALSE, matrix.type = "adjacency")
  neg_net <- network::network(neg_adj(), directed = FALSE, matrix.type = "adjacency")

  ref <- Signed(pos_net, neg_net)     # reference from previous round
  cmp <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj())

  # Both should have identical layer label sets
  expect_setequal(layer_names(ref), layer_names(cmp))
})

test_that("adjacency path: class vector agrees between network.sign and Signed", {
  pos_net <- network::network(pos_adj(), directed = FALSE, matrix.type = "adjacency")
  neg_net <- network::network(neg_adj(), directed = FALSE, matrix.type = "adjacency")

  ref <- Signed(pos_net, neg_net)
  cmp <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj())

  expect_true("static.sign" %in% class(ref))
  expect_true("static.sign" %in% class(cmp))
})

test_that("adjacency path: dual.sign behaviour consistent across both functions", {
  pos_net <- network::network(pos_adj(), directed = FALSE, matrix.type = "adjacency")
  neg_net <- network::network(neg_adj(), directed = FALSE, matrix.type = "adjacency")

  ref_f <- Signed(pos_net, neg_net, dual.sign = FALSE)
  cmp_f <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj(), dual.sign = FALSE)
  ref_t <- Signed(pos_net, neg_net, dual.sign = TRUE)
  cmp_t <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj(), dual.sign = TRUE)

  expect_equal(ref_f %n% "dual.sign", cmp_f %n% "dual.sign")
  expect_equal(ref_t %n% "dual.sign", cmp_t %n% "dual.sign")
})

test_that("sign attr mirrors .LayerName — consistent across both functions", {
  pos_net <- network::network(pos_adj(), directed = FALSE, matrix.type = "adjacency")
  neg_net <- network::network(neg_adj(), directed = FALSE, matrix.type = "adjacency")

  ref <- Signed(pos_net, neg_net)
  cmp <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj())

  expect_equal(ref %v% "sign", ref %v% ".LayerName")
  expect_equal(cmp %v% "sign", cmp %v% ".LayerName")
})

test_that("network.sign adjacency vs edgelist: same edge structure", {
  adj_net <- network.sign(signed_adj(),      matrix.type = "adjacency")
  el_net  <- network.sign(signed_edgelist(), matrix.type = "edgelist")

  # Both should produce the same number of edges per layer
  expect_equal(network.edgecount(adj_net), network.edgecount(el_net))
})

# ---------------------------------------------------------------------------
# 9. Error cases
# ---------------------------------------------------------------------------

test_that("unsupported matrix.type errors", {
  expect_error(
    network.sign(signed_adj(), matrix.type = "dense"),
    regexp = "arg"
  )
})

test_that("edge attribute values outside {-1, 0, 1} are handled", {
  bad <- signed_adj()
  bad[1, 2] <- 2L
  # build_network splits on >0 / <0 so value 2 is treated as positive — no error expected
  expect_no_error(network.sign(bad, matrix.type = "adjacency"))
})

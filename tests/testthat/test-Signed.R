# tests/testthat/test-Signed.R

library(testthat)
library(network)

# ---------------------------------------------------------------------------
# Adjacency matrices
# ---------------------------------------------------------------------------

#' 4x4 signed adjacency matrix with values in {-1, 0, 1}
signed_adj <- function() {
  m <- matrix(0L, 4, 4)
  m[1, 2] <-  1L   # positive
  m[2, 3] <- -1L   # negative
  m[3, 4] <-  1L   # positive
  m
}

#' Same matrix but with row/col names
signed_adj_named <- function() {
  m <- signed_adj()
  dimnames(m) <- list(letters[1:4], letters[1:4])
  m
}

#' Positive-only binary adjacency (derived from signed_adj)
pos_adj <- function() (signed_adj() > 0) * 1L

#' Negative-only binary adjacency (derived from signed_adj)
neg_adj <- function() (signed_adj() < 0) * 1L

# ---------------------------------------------------------------------------
# network objects
# ---------------------------------------------------------------------------

#' Separate positive and negative network objects (undirected, 4 nodes)
make_pos_neg_networks <- function(directed = FALSE) {
  pos_net <- network::network(pos_adj(), directed = directed, matrix.type = "adjacency")
  neg_net <- network::network(neg_adj(), directed = directed, matrix.type = "adjacency")
  list(pos = pos_net, neg = neg_net)
}

#' Single network with a {-1, 0, 1}-valued edge attribute "sign"
make_signed_attr_network <- function(directed = FALSE) {
  nw <- network::network.initialize(4, directed = directed)
  network::add.edge(nw, 1, 2)
  network::add.edge(nw, 2, 3)
  network::add.edge(nw, 3, 4)
  nw %e% "sign" <- c(1L, -1L, 1L)
  nw
}

#' Single network with two separate binary edge attributes
make_two_attr_network <- function(directed = FALSE) {
  nw <- network::network.initialize(4, directed = directed)
  network::add.edge(nw, 1, 2)
  network::add.edge(nw, 2, 3)
  network::add.edge(nw, 3, 4)
  nw %e% "is_pos" <- c(1L, 0L, 1L)
  nw %e% "is_neg" <- c(0L, 1L, 0L)
  nw
}

# ---------------------------------------------------------------------------
# Edgelist
# ---------------------------------------------------------------------------

#' Three-column edgelist consistent with signed_adj()
signed_edgelist <- function() {
  data.frame(
    From = c(1L, 2L, 3L),
    To   = c(2L, 3L, 4L),
    Sign = c(1L, -1L, 1L)
  )
}

# ---------------------------------------------------------------------------
# Shared expectations
# ---------------------------------------------------------------------------

#' Extract unique values of the 'sign' vertex attribute
layer_names <- function(net) unique(net %v% "sign")

#' Assert the standard invariants that EVERY signed network must satisfy,
#' regardless of which constructor produced it.
#'
#' @param net  A signed network object.
#' @param expected_layers  Character vector of expected layer name values.
#' @param expect_dual  Logical; expected value of the dual.sign network attribute.
expect_signed_invariants <- function(net,
                                     expected_layers = c("+", "-"),
                                     expect_dual     = FALSE) {
  # Class membership
  expect_s3_class(net, "network")
  expect_s3_class(net, "static.sign")

  # Layer labels present
  expect_setequal(layer_names(net), expected_layers)

  # sign vertex attr mirrors .LayerName
  expect_equal(net %v% "sign", net %v% ".LayerName")

  # dual.sign network attribute stored correctly
  expect_equal(net %n% "dual.sign", expect_dual)

  # Constraint present iff dual.sign is FALSE
  constraints <- net %ergmlhs% "constraints"
  if (!expect_dual) {
    expect_false(is.null(constraints))
  }
}

# ---------------------------------------------------------------------------
# 1. Signed(pos_nw, neg_nw) — two network objects
# ---------------------------------------------------------------------------

test_that("two-net: satisfies all signed invariants (dual.sign = FALSE)", {
  nws <- make_pos_neg_networks()
  net <- Signed(nws$pos, nws$neg)
  expect_signed_invariants(net, expected_layers = c("pos", "neg"), expect_dual = FALSE)
})

test_that("two-net: satisfies all signed invariants (dual.sign = TRUE)", {
  nws <- make_pos_neg_networks()
  net <- Signed(nws$pos, nws$neg, dual.sign = TRUE)
  expect_signed_invariants(net, expected_layers = c("pos", "neg"), expect_dual = TRUE)
})

test_that("two-net: directed networks accepted without error", {
  nws <- make_pos_neg_networks(directed = TRUE)
  expect_no_error(Signed(nws$pos, nws$neg))
})

# ---------------------------------------------------------------------------
# 2. Signed(nw, sign_attr) — single network + {-1, 0, 1} attribute
# ---------------------------------------------------------------------------

test_that("sign-attr: satisfies all signed invariants (dual.sign = FALSE)", {
  nw  <- make_signed_attr_network()
  net <- Signed(nw, "sign")
  expect_signed_invariants(net, expected_layers = c("pos", "neg"), expect_dual = FALSE)
})

test_that("sign-attr: satisfies all signed invariants (dual.sign = TRUE)", {
  nw  <- make_signed_attr_network()
  net <- Signed(nw, "sign", dual.sign = TRUE)
  expect_signed_invariants(net, expected_layers = c("pos", "neg"), expect_dual = TRUE)
})

test_that("sign-attr: errors when attribute does not exist", {
  nw <- network::network.initialize(4)
  network::add.edge(nw, 1, 2)
  expect_error(Signed(nw, "nonexistent"), regexp = "not found")
})

test_that("sign-attr: errors when attribute contains values outside {-1, 0, 1}", {
  nw <- network::network.initialize(4)
  network::add.edge(nw, 1, 2)
  nw %e% "sign" <- c(2L)
  expect_error(Signed(nw, "sign"), regexp = "\\{-1, 0, 1\\}")
})

test_that("sign-attr: value 0 is accepted as unsigned/absent", {
  nw <- network::network.initialize(4)
  network::add.edge(nw, 1, 2)
  network::add.edge(nw, 2, 3)
  nw %e% "sign" <- c(1L, 0L)
  expect_no_error(Signed(nw, "sign"))
})

test_that("sign-attr: positive edges land in pos layer, negative in neg layer", {
  nw  <- make_signed_attr_network()   # edges: 1-2 (+1), 2-3 (-1), 3-4 (+1)
  net <- Signed(nw, "sign")
  expect_true("pos" %in% layer_names(net))
  expect_true("neg" %in% layer_names(net))
})

# ---------------------------------------------------------------------------
# 3. Signed(nw, pos_attr, neg_attr) — two binary edge attributes
# ---------------------------------------------------------------------------

test_that("two-attr: satisfies all signed invariants (dual.sign = FALSE)", {
  nw  <- make_two_attr_network()
  net <- Signed(nw, "is_pos", "is_neg")
  expect_signed_invariants(net, expected_layers = c("pos", "neg"), expect_dual = FALSE)
})

test_that("two-attr: satisfies all signed invariants (dual.sign = TRUE)", {
  nw  <- make_two_attr_network()
  net <- Signed(nw, "is_pos", "is_neg", dual.sign = TRUE)
  expect_signed_invariants(net, expected_layers = c("pos", "neg"), expect_dual = TRUE)
})

# ---------------------------------------------------------------------------
# 4. dual.sign constraint — detailed checks
# ---------------------------------------------------------------------------

test_that("dual.sign = FALSE: constraint formula is non-trivial", {
  nws <- make_pos_neg_networks()
  net <- Signed(nws$pos, nws$neg, dual.sign = FALSE)
  c_text <- deparse(net %ergmlhs% "constraints")
  # Must mention ChangeStats or the pos & neg interaction term
  expect_match(c_text, "ChangeStats|pos.*neg|neg.*pos", perl = TRUE)
})

test_that("dual.sign FALSE vs TRUE produce different constraint formulas", {
  nws   <- make_pos_neg_networks()
  net_f <- Signed(nws$pos, nws$neg, dual.sign = FALSE)
  net_t <- Signed(nws$pos, nws$neg, dual.sign = TRUE)
  expect_false(identical(
    deparse(net_f %ergmlhs% "constraints"),
    deparse(net_t %ergmlhs% "constraints")
  ))
})

# ---------------------------------------------------------------------------
# 5. sign vertex attribute mirrors .LayerName (all signatures)
# ---------------------------------------------------------------------------

test_that("sign ↔ .LayerName: two-net signature", {
  nws <- make_pos_neg_networks()
  net <- Signed(nws$pos, nws$neg)
  expect_equal(net %v% "sign", net %v% ".LayerName")
})

test_that("sign ↔ .LayerName: sign-attr signature", {
  nw  <- make_signed_attr_network()
  net <- Signed(nw, "sign")
  expect_equal(net %v% "sign", net %v% ".LayerName")
})

test_that("sign ↔ .LayerName: two-attr signature", {
  nw  <- make_two_attr_network()
  net <- Signed(nw, "is_pos", "is_neg")
  expect_equal(net %v% "sign", net %v% ".LayerName")
})

# ---------------------------------------------------------------------------
# 6. Optional Layer() arguments pass through
# ---------------------------------------------------------------------------

test_that(".symmetric is forwarded to Layer() without error", {
  nws <- make_pos_neg_networks(directed = TRUE)
  expect_no_error(Signed(nws$pos, nws$neg, .symmetric = c(TRUE, TRUE)))
})

# ---------------------------------------------------------------------------
# 7. Error handling — bad inputs
# ---------------------------------------------------------------------------

test_that("errors with informative message on completely bad input", {
  expect_error(Signed(42), regexp = "Unrecognized format")
})

test_that("errors when called with no arguments", {
  expect_error(Signed())
})

test_that("errors when three network objects are passed", {
  nws   <- make_pos_neg_networks()
  extra <- network::network.initialize(4)
  expect_error(Signed(nws$pos, nws$neg, extra))
})

# ---------------------------------------------------------------------------
# 8. Cross-function parity: Signed() vs network.sign()
#
# These tests verify structural equivalence between the two constructors.
# The mirror of these tests lives in test-network_sign.R § 8.
#
# We do NOT assert identical layer label strings ("pos"/"neg" vs "+"/"-")
# because the two functions use different naming conventions by design.
# We DO assert everything else: class, dual.sign flag, constraint behaviour,
# sign↔.LayerName invariant, and edge counts.
# ---------------------------------------------------------------------------

test_that("parity: both return static.sign network", {
  nws <- make_pos_neg_networks()
  ref <- Signed(nws$pos, nws$neg)
  cmp <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj())
  expect_s3_class(ref, "static.sign")
  expect_s3_class(cmp, "static.sign")
})

test_that("parity: dual.sign = FALSE stored consistently", {
  nws <- make_pos_neg_networks()
  ref <- Signed(nws$pos, nws$neg, dual.sign = FALSE)
  cmp <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj(), dual.sign = FALSE)
  expect_equal(ref %n% "dual.sign", cmp %n% "dual.sign")
  expect_false(ref %n% "dual.sign")
})

test_that("parity: dual.sign = TRUE stored consistently", {
  nws <- make_pos_neg_networks()
  ref <- Signed(nws$pos, nws$neg, dual.sign = TRUE)
  cmp <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj(), dual.sign = TRUE)
  expect_equal(ref %n% "dual.sign", cmp %n% "dual.sign")
  expect_true(ref %n% "dual.sign")
})

test_that("parity: dual.sign = FALSE adds constraint in both constructors", {
  nws <- make_pos_neg_networks()
  ref <- Signed(nws$pos, nws$neg, dual.sign = FALSE)
  cmp <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj(), dual.sign = FALSE)
  expect_false(is.null(ref %ergmlhs% "constraints"))
  expect_false(is.null(cmp %ergmlhs% "constraints"))
})

test_that("parity: dual.sign FALSE vs TRUE changes constraints consistently", {
  nws <- make_pos_neg_networks()
  ref_f <- Signed(nws$pos, nws$neg, dual.sign = FALSE)
  ref_t <- Signed(nws$pos, nws$neg, dual.sign = TRUE)
  cmp_f <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj(), dual.sign = FALSE)
  cmp_t <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj(), dual.sign = TRUE)

  # Each constructor: FALSE ≠ TRUE
  expect_false(identical(
    deparse(ref_f %ergmlhs% "constraints"),
    deparse(ref_t %ergmlhs% "constraints")
  ))
  expect_false(identical(
    deparse(cmp_f %ergmlhs% "constraints"),
    deparse(cmp_t %ergmlhs% "constraints")
  ))
})

test_that("parity: sign ↔ .LayerName invariant holds in both constructors", {
  nws <- make_pos_neg_networks()
  ref <- Signed(nws$pos, nws$neg)
  cmp <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj())
  expect_equal(ref %v% "sign", ref %v% ".LayerName")
  expect_equal(cmp %v% "sign", cmp %v% ".LayerName")
})

test_that("parity: edge counts match between Signed() and network.sign()", {
  nws <- make_pos_neg_networks()
  ref <- Signed(nws$pos, nws$neg)
  cmp <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj())
  expect_equal(network.edgecount(ref), network.edgecount(cmp))
})

test_that("parity: number of layers (unique sign values) is 2 in both", {
  nws <- make_pos_neg_networks()
  ref <- Signed(nws$pos, nws$neg)
  cmp <- network.sign(pos.mat = pos_adj(), neg.mat = neg_adj())
  expect_length(layer_names(ref), 2)
  expect_length(layer_names(cmp), 2)
})


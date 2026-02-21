library(testthat)
library(ergm)

testthat::test_that("timepoint pooling works correctly (majority + tie breaker)", {

  # helper: build adjacency matrices
  make_mat <- function(edges, n = 3) {
    m <- matrix(0, n, n)
    for (e in edges) {
      m[e[1], e[2]] <- e[3]
    }
    m
  }

  # Create 6 timepoint matrices
  # We pool into: (1:3) and (4:6)
  mats <- list(
    make_mat(list(c(1,2, 1), c(2,3,-1))),  # t1
    make_mat(list(c(1,2, 1), c(2,3,-1))),  # t2
    make_mat(list(c(1,2,-1), c(2,3,-1))),  # t3  -> makes (1,2) tie (2 pos vs 1 neg)
    make_mat(list(c(1,2,-1), c(2,3, 1))),  # t4
    make_mat(list(c(1,2,-1), c(2,3, 1))),  # t5
    make_mat(list(c(1,2,-1), c(2,3,-1)))   # t6
  )

  # Define pooling windows
  timepoints <- list(1:3, 4:6)

  # Run pooling through network.sign()
  # tie.breaker = "zero" means ties become 0
  pooled <- network.sign(
    mat = mats,
    matrix.type = "adjacency",
    timepoints = timepoints,
    tie.breaker = "zero"
  )

  testthat::expect_equal(length(pooled$gal$NetList), 2)

  # Check names exist and match expected pooling windows
  testthat::expect_equal(pooled$gal$names, c("Timepoints 1-3", "Timepoints 4-6"))

  # Convert back to adjacency to check pooling result
  get_signed_adj <- function(nets) {
    tmp <- UnLayer(nets)
    lapply(tmp, function(net) as.sociomatrix(net, attrname = "sign"))
  }
  nets <- UnLayer(pooled)
  pooled1 <- as.sociomatrix(nets[[1]], attrname = "sign")
  pooled2 <- as.sociomatrix(nets[[2]], attrname = "sign")

  # --- Check pooled window 1: timepoints 1:3 ---
  # Edge (1,2): +1, +1, -1  => pos=2 neg=1 => +1
  testthat::expect_equal(pooled1[1,2], 1)

  # Edge (2,3): -1, -1, -1 => neg dominates => -1
  testthat::expect_equal(pooled1[2,3], -1)

  # --- Check pooled window 2: timepoints 4:6 ---
  # Edge (1,2): -1, -1, -1 => -1
  testthat::expect_equal(pooled2[1,2], -1)

  # Edge (2,3): +1, +1, -1 => +1
  testthat::expect_equal(pooled2[2,3], 1)

})

testthat::test_that("tie breaker works correctly", {
  suppressWarnings({
  make_mat <- function(edges, n = 2) {
    m <- matrix(0, n, n)
    for (e in edges) {
      m[e[1], e[2]] <- e[3]
    }
    m
  }

  # Construct a real tie:
  # Edge (1,2) is +1 in t1 and -1 in t2 -> pos=1 neg=1 => tie
  mats <- list(
    make_mat(list(c(1,2, 1))),   # t1
    make_mat(list(c(1,2,-1)))    # t2
  )

  timepoints <- 1

  # --------------------------
  # tie.breaker = "zero"
  pooled_zero <- network.sign(
    mat = mats,
    matrix.type = "adjacency",
    timepoints = timepoints,
    tie.breaker = "zero"
  )

  nets_zero <- UnLayer(pooled_zero)
  pooled_mat_zero <- as.sociomatrix(nets_zero, attrname = "sign")

  testthat::expect_equal(pooled_mat_zero[1,2], 0)


  # --------------------------
  # tie.breaker = "positive"
  pooled_pos <- network.sign(
    mat = mats,
    matrix.type = "adjacency",
    timepoints = timepoints,
    tie.breaker = "positive"
  )

  nets_pos <- UnLayer(pooled_pos)
  pooled_mat_pos <- as.sociomatrix(nets_pos, attrname = "sign")

  testthat::expect_equal(pooled_mat_pos[1,2], 1)


  # --------------------------
  # tie.breaker = "negative"
  pooled_neg <- network.sign(
    mat = mats,
    matrix.type = "adjacency",
    timepoints = timepoints,
    tie.breaker = "negative"
  )

  nets_neg <- UnLayer(pooled_neg)
  pooled_mat_neg <- as.sociomatrix(nets_neg, attrname = "sign")

  testthat::expect_equal(pooled_mat_neg[1,2], -1)


  # --------------------------
  # tie.breaker = "first"
  pooled_first <- network.sign(
    mat = mats,
    matrix.type = "adjacency",
    timepoints = timepoints,
    tie.breaker = "first"
  )

  nets_first <- UnLayer(pooled_first)
  pooled_mat_first <- as.sociomatrix(nets_first, attrname = "sign")

  # first timepoint was +1
  testthat::expect_equal(pooled_mat_first[1,2], 1)


  # --------------------------
  # tie.breaker = "last"
  pooled_last <- network.sign(
    mat = mats,
    matrix.type = "adjacency",
    timepoints = timepoints,
    tie.breaker = "last"
  )

  nets_last <- UnLayer(pooled_last)
  pooled_mat_last <- as.sociomatrix(nets_last, attrname = "sign")

  # last timepoint was -1
  testthat::expect_equal(pooled_mat_last[1,2], -1)
  })
})

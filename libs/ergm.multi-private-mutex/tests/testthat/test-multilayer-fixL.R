#  File tests/testthat/test-multilayer-MLE.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
nsim <- 5
nw0 <- network.initialize(5, dir=FALSE)

test_net <- function(nw0, fix, max_invariant, invariant){
  n <- network.size(nw0)
  d <- network.dyadcount(nw0)
  nwl0 <- Layer(nw0, nw0)

  # Maximum test
  nwl1 <- simulate(nwl0 ~ L(~edges, ~`1`) + L(~edges, ~`2`), constraints=~.+fixL(fix), coef=c(Inf,Inf), nsim=nsim)

  expect_equal(summary(nwl1~L(~edges, fix)), numeric(nsim), ignore_attr=TRUE)
  expect_equal(cbind(d, summary(nwl1~L(~edges, ~`1`) + L(~edges, ~`2`) + L(~edges, ~`1`&`2`)))%*%max_invariant, numeric(nsim), ignore_attr=TRUE)

  # random test
  nwl1 <- simulate(nwl0 ~ L(~edges, ~`1`) + L(~edges, ~`2`), constraints=~.+fixL(fix), coef=c(0,0), nsim=nsim)

  expect_equal(summary(nwl1~L(~edges, fix)), numeric(nsim), ignore_attr=TRUE)
  expect_equal(summary(nwl1~L(~edges, ~`1`) + L(~edges, ~`2`) + L(~edges, ~`1`&`2`))%*%invariant, numeric(nsim), ignore_attr=TRUE)
}

test_net(nw0, ~`1`&`2`, c(-1,1,1,0), c(0,0,1))
test_net(nw0, ~`1`&!`2`, c(-2,1,1,0), c(-1,0,1))

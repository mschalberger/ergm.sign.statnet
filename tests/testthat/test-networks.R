test_that("multinetwork network representation works", {
  set.seed(1)

  adj_matrix_1 <- matrix(sample(c(-1, 0, 1), 100, replace = TRUE), nrow = 10)
  adj_matrix_2 <- matrix(sample(c(-1, 0, 1), 25, replace = TRUE), nrow = 5)

  full_adj <- rbind(cbind(adj_matrix_1, matrix(0, nrow = nrow(adj_matrix_1), ncol = ncol(adj_matrix_2))),
        cbind(matrix(0, nrow = nrow(adj_matrix_2), ncol = ncol(adj_matrix_1)), adj_matrix_2))

  full_net <- network.sign(full_adj, matrix.type = "adjacency", directed = T)
  set.vertex.attribute(full_net, "adj", c(rep(1, nrow(adj_matrix_1)), rep(2, nrow(adj_matrix_2))))

  net1 <- network.sign(adj_matrix_1, matrix.type = "adjacency", directed = T)
  net2 <- network.sign(adj_matrix_2, matrix.type = "adjacency", directed = T)

  # fit ergm
  s1 <- summary(full_net ~ Pos(~edges) + Neg(~edges))
  s2 <- summary(Networks(full_net) ~N(~Pos(~edges) + Neg(~edges)))

  expect_equal(unname(s1),unname(s2))
})

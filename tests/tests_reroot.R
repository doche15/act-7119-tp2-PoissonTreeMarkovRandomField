AA <- matrix(c(1, "alpha12", "alpha13", 0, 0, 0, 0,
               "alpha12", 1, 0, 0, 0, 0, 0,
               "alpha13", 0, 1, "alpha34", "alpha35", 0, 0,
               0, 0, "alpha34", 1, 0, "alpha46", "alpha47",
               0, 0, "alpha35", 0, 1, 0, 0,
               0, 0, 0, "alpha46", 0, 1, 0,
               0, 0, 0, "alpha47", 0, 0, 1), byrow = TRUE, ncol = 7)



reroot(AA, 1, 3)

alpha12 <- 0.1
alpha13 <- 0.2
alpha34 <- 0.4
alpha35 <- 0.6
alpha46 <- 0.3
alpha47 <- 0.7

AA <- matrix(c(1, alpha12, alpha13, 0, 0, 0, 0,
               alpha12, 1, 0, 0, 0, 0, 0,
               alpha13, 0, 1, alpha34, alpha35, 0, 0,
               0, 0, alpha34, 1, 0, alpha46, alpha47,
               0, 0, alpha35, 0, 1, 0, 0,
               0, 0, 0, alpha46, 0, 1, 0,
               0, 0, 0, alpha47, 0, 0, 1), byrow = TRUE, ncol = 7)
lam <- 50
rea <- rMPMRF(10000000, AA, lam)
reaM <- rowSums(rea)

mean(rea[, 7] * (reaM == 400))

exp_alloc_Nv(AA, 1, lam, 7)[401]


reorder_adjacency <- function(A, r, r_prime) {
    d <- nrow(A)

    for (k in seq(d, 2, by = -1)) {
        pi_k <- which(A[k, ] > 0)[1]  # Find the first nonzero column in row k
    }

    w <- r_prime
    w_vector <- c()
    while (w != r) {
        w_vector <- c(w_vector, w)
        w <- pi_k  # Overwrite w by pi_k
    }

    for (i in w_vector) {
        for (j in seq(i - 1, 1, by = -1)) {
            A[c(j, j + 1), ] <- A[c(j + 1, j), ]  # Switch rows j and j+1
            A[, c(j, j + 1)] <- A[, c(j + 1, j)]  # Switch columns j and j+1
        }
    }

    return(A)
}

reorient_tree_igraph <- function(adj_matrix, new_root) {
    # Create an undirected graph from the adjacency matrix
    d <- nrow(adj_matrix)
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected",
                                     weighted = TRUE, diag = FALSE)

    # Perform BFS from the new root to determine the new direction of edges
    bfs_result <- bfs(g, root = new_root, neimode = "all", order = TRUE, father = TRUE)

    # Extract the parent-child relationships
    parent <- as.numeric(bfs_result$father)
    parent[is.na(parent)] <- new_root  # Ensure root has itself as a parent

    # Create a new adjacency matrix with directed edges
    new_adj_matrix <- diag(d)
    for (i in 1:nrow(adj_matrix)) {
        for (j in 1:nrow(adj_matrix)) {
            new_adj_matrix[i, j] <- adj_matrix[bfs_result$order[i], bfs_result$order[j]]# weight
        }
    }

    return(new_adj_matrix)
}

reorient_tree_igraph(A, 3)

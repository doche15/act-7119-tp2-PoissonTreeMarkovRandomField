###
### ACT-7119, travail 2
### Tests JG, algorithme 5
###
##


changer_racine_arbre <- function(A, nouvelle_racine, ancienne_racine = 1)
{
    dimension <- nrow(A)

    # m <- 0
    # pi_k <- numeric(dimension - 1)
    # for (i in rev(2:d))
    # {
    #     m <- m + 1
    #     pi_k[m] <- min(which(A[i,] > 0))
    # }

    ww <- numeric()
    w <- nouvelle_racine
    while (w != ancienne_racine) {
        print(w)
        ww <- c(ww, w)
        w <- min(which(A[w,] > 0))
    }
    print(ww)

    nouveau_A <- A
    for (i in ww)
    {
        for (j in rev(1:(i - 1)))
        {
            # Changements, lignes
            temp_lignej <- nouveau_A[j,]

            nouveau_A[j,] <- nouveau_A[j + 1,]
            nouveau_A[j + 1,] <- temp_lignej

            # Changements, colonnes
            temp_colonnej <- nouveau_A[, j]

            nouveau_A[, j] <- nouveau_A[, j + 1]
            nouveau_A[, j + 1] <- temp_colonnej
        }
    }

    nouveau_A
}


reroot <- function(adj_matrix, new_root) {

    d <- nrow(adj_matrix)
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected",
                                     weighted = TRUE, diag = FALSE)

    # Chemin de la racine aux autres sommets
    bfs_result <- bfs(g, root = new_root, neimode = "all", order = TRUE, father = TRUE)
    print(bfs_result)


    # parent <- as.numeric(bfs_result$father)
    # parent[is.na(parent)] <- new_root

    # Élément par élément créer la matrice
    new_adj_matrix <- diag(d)
    for (i in 1:nrow(adj_matrix)) {
        for (j in 1:nrow(adj_matrix)) {
            new_adj_matrix[i, j] <- adj_matrix[bfs_result$order[i], bfs_result$order[j]]# weight
        }
    }

    return(new_adj_matrix)
}


a12 <- 0.2
a13 <- 0.5
a34 <- 0.3
a35 <- 0.47
a46 <- 0.33
a47 <- 0.4
A <- matrix(c(1,   a12, a13, 0,   0,   0,   0,
              a12, 1,   0,   0,   0,   0,   0,
              a13, 0,   1,   a34, a35, 0,   0,
              0,   0,   a34, 1,   0,   a46, a47,
              0,   0,   a35, 0,   1,   0,   0,
              0,   0,   0,   a46, 0,   1,   0,
              0,   0,   0,   a47, 0,   0,   1),
            nrow = 7,
            byrow = TRUE)


changer_racine_arbre(A, 3)
reroot(A, 3)


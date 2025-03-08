###
### Travail 2, ACT-7119
### Algorithm 5 : Changing root in accordance to with an adjacency matrix is in topological order.

reroot = function(A, root_node, new_root){
    # A : matrice d'adjacence
    # root_node : racine initiale
    # new_root : nouvelle racine voulue

    d = nrow(A)

    omega = new_root

    omega_vec = numeric()
    parents <- numeric(d - 1)
    for (k in d:2)
        parents[k - 1] <- min(which(A[k, ] > 0))


    while (omega != root_node){

        omega_vec = c(omega_vec, omega)

        omega = parents[omega - 1]

    }
    omega_vec <- sort(omega_vec)
    for (i in omega_vec){

        for (j in rev(seq(i - 1))){

            row_j = A[j, ] # ligne j
            row_jj = A[j + 1, ] # ligne j + 1

            # Switch
            A[j, ] = row_jj
            A[j + 1, ] = row_j

            col_j = A[, j] # colonne j
            col_jj = A[, j + 1] # colonne j + 1

            # Switch
            A[, j] = col_jj
            A[, j + 1] = col_j

        }

    }

    A

}
#
# AA <- matrix(c(1, "alpha12", "alpha13", 0, 0, 0, 0,
#                "alpha12", 1, 0, 0, 0, 0, 0,
#                "alpha13", 0, 1, "alpha34", "alpha35", 0, 0,
#                0, 0, "alpha34", 1, 0, "alpha46", "alpha47",
#                0, 0, "alpha35", 0, 1, 0, 0,
#                0, 0, 0, "alpha46", 0, 1, 0,
#                0, 0, 0, "alpha47", 0, 0, 1), byrow = TRUE, ncol = 7)
#
# AA_prime <- reroot(AA, 1, 3)
#
#
# alpha12 <- 0.1
# alpha13 <- 0.2
# alpha34 <- 0.4
# alpha35 <- 0.6
# alpha46 <- 0.3
# alpha47 <- 0.7
#
# AA <- matrix(c(1, alpha12, alpha13, 0, 0, 0, 0,
#                alpha12, 1, 0, 0, 0, 0, 0,
#                alpha13, 0, 1, alpha34, alpha35, 0, 0,
#                0, 0, alpha34, 1, 0, alpha46, alpha47,
#                0, 0, alpha35, 0, 1, 0, 0,
#                0, 0, 0, alpha46, 0, 1, 0,
#                0, 0, 0, alpha47, 0, 0, 1), byrow = TRUE, ncol = 7)
#
# AA_pr <- reroot(AA, 1, 3)
#
# reaA <- rMPMRF(10000000, AA, lam)
# reaApr <- rMPMRF(10000000, AA_pr, lam)
#
# colMeans(reaA)
# colMeans(reaApr)
#
#
# reaMA <- rowSums(reaA)
# reaMApr <- rowSums(reaApr)
#
# mean(reaMA)
# mean(reaMApr)
#
# mean(reaMA^2) - mean(reaMA)^2
# mean(reaMApr^2) - mean(reaMApr)^2
#
#
# ccc <- exp_alloc_Nv(AA, 1, lam, 4)
#
# sum(ccc) # somme des allocations donne bien l'espérance
#
#
# mean(reaA[, 4] * (reaMA == 10))
# mean(reaApr[, 1] * (reaMApr == 10)) # avec matrice ajustée
# ccc[11]
#
# cccAppr <- exp_alloc_Nv(AA_pr, 1, lam, 4)
# ccc[11] ## Semble ok

## Attention: les résultats étaient moins proches avec mon arbre de d=32,
##            mais peut-être que la variance est trop élevée, ce qui fait que ça diverge un peu


###
### Travail 2, ACT-7119
### Exemple
###
library(igraph)
d <- 32
A <- diag(32) # Matrice identité de d = 32

## Arètes
# (1, 2), (2, 3), (2, 4), (2, 5), (5, 6), (5, 7),
# (6, 8), (7, 9), (9, 11), (9, 14), (9, 13), (9, 12),
# (9, 11), (9, 10), (1, 16), (16, 18), (16, 17), (16, 19),
# (16, 20), (20, 22), (22, 21), (21, 23), (23, 24), (23, 25),
# (23, 26), (22, 28), (28, 27), (28, 29), (28, 30), (28, 31), (31, 32)

aretes <- matrix(c(1, 2, 2, 3, 2, 4, 2, 5, 5, 6, 5, 7,
                   6, 8, 7, 9, 9, 11, 9, 14, 9, 13, 9, 12,
                   9, 11, 9, 10, 1, 16, 16, 18, 16, 17, 16, 19,
                   16, 20, 20, 22, 22, 21, 21, 23, 23, 24, 23, 25,
                   23, 26, 22, 28, 28, 27, 28, 29, 28, 30, 28, 31, 31, 32),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, c(0.2, 0.1, 0.1, 0.7, 0.3, 0.3, 0.4, 0.4, rep(0.6, 6),
                 0.7, 0.02, 0.03, 0.04, 0.5, 0.5, 0.4, 0.4, rep(0.1, 3),
                 0.9, rep(0.7, 4), 0.95))

for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A[u, v] <- aretes_alpha[i, 3]
    A[v, u] <- aretes_alpha[i, 3]
}

adj <- A
lam <- 4

source("graph_A.R")
see_tree_graph(adj, root_node = 1, titre = "Arbre d = 32")


## Simulation
source("sample_MPMRF.R")
nn <- 10000000
realisations <- rMPMRF(nn, adj, lam)

MM <- rowSums(realisations)


mean(MM)
mean(MM^2) - mean(MM)^2

adjpr <- reroot(adj, 1, 4)

realisationspr <- rMPMRF(nn, adjpr, lam)

MMpr <- rowSums(realisationspr)


mean(MMpr)
mean(MMpr^2) - mean(MMpr)^2

exp_alloc_Nv(adj, 1, lam, 4)[100]
exp_alloc_Nv(adjpr, 1, lam, 1)[100]

mean(realisations[, 4] * (MM == 99))

mean(realisationspr[, 1] * (MMpr == 99))

## Exact
source("pdf_M.R")

nfft <- 2^9
fM <- pdf_M(adj, lam, nfft)
sum(fM)
sum((0:(nfft - 1)) * fM)
sum((0:(nfft - 1))^2 * fM) - sum((0:(nfft - 1)) * fM)^2

cumsum(fM)[90 + 1]
mean(MM <= 90)

cumsum(fM)[98 + 1]
mean(MM <= 98)

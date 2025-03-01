library(igraph)
source("graph_A.R")
source("sample_MPMRF.R")
source("pdf_M.R")

alpha12 <- 0.4
alpha13 <- 0.6
alpha34 <- 0.5
alpha35 <- 0.3
alpha46 <- 0.7
alpha47 <- 0.8
lam <- 6
adj <- matrix(c(1, alpha12, alpha13, 0, 0, 0, 0,
                alpha12, 1, 0, 0, 0, 0, 0,
                alpha13, 0, 1, alpha34, alpha35, 0, 0,
                0, 0, alpha34, 1, 0, alpha46, alpha47,
                0, 0, alpha35, 0, 1, 0, 0,
                0, 0, 0, alpha46, 0, 1, 0,
                0, 0, 0, alpha47, 0, 0, 1), ncol = 7, byrow = TRUE)

see_tree_graph(adj, 1, titre = "Figure 2")

nn <- 300000000
realisations <- rMPMRF(nn, adj, lam)
reaN <- rowSums(realisations)
mean(reaN)
mean(reaN^2) - mean(reaN)^2

colMeans(realisations)

apply(realisations, 2, var)
cor(realisations[, 1], realisations[, 2])
cor(realisations[, 2], realisations[, 3])
cor(realisations[, 2], realisations[, 4])
cor(realisations[, 1], realisations[, 3])
cor(realisations[, 2], realisations[, 7])

nfft <- 2^15

fM <- pdf_M(adj, lam)
sum(fM)

sum((0:(nfft - 1)) * fM)
sum((0:(nfft - 1))^2 * fM) - sum((0:(nfft - 1)) * fM)^2

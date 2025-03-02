###
### ACT-7119
### Test arbe d = 4
### [Côté et al., 2025].
###
##

library(igraph)
source("graph_A.R")

alpha12 <- 0.4
alpha23 <- 0.5
alpha24 <- 0.6
adj <- matrix(c(1, alpha12, 0, 0,
              alpha12, 1, alpha23, alpha24,
              0, alpha23, 1, 0,
              0, alpha24, 0, 1), ncol = 4, byrow = TRUE)

see_tree_graph(adj, root_node = 1, titre = "Graphique d = 4")

lam <- 5

## Simulation à la main

rea <- function(a12, a23, a24, lambda)
{
    N1 <- rpois(1, lambda)
    prop2 <- rbinom(1, N1, a12)
    inno2 <- rpois(1, lambda * (1 - a12))

    N2 <- prop2 + inno2

    prop3 <- rbinom(1, N2, a23)
    inno3 <- rpois(1, lambda * (1 - a23))

    N3 <- prop3 + inno3

    prop4 <- rbinom(1, N2, a24)
    inno4 <- rpois(1, lambda * (1 - a24))

    N4 <- prop4 + inno4

    c(N1, N2, N3, N4)
}

set.seed(20250301)
nn <- 100000
realisations <- t(replicate(nn, rea(alpha12, alpha23, alpha24, lam)))
reaN <- rowSums(realisations)

colMeans(realisations)
apply(realisations, 2, var)
cor(realisations[, 1], realisations[, 2])
cor(realisations[, 2], realisations[, 3])
cor(realisations[, 2], realisations[, 4])
cor(realisations[, 1], realisations[, 3])
cor(realisations[, 1], realisations[, 4])

4 * lam + 2 * alpha12 * lam + 2 * alpha23 * lam + 2 * alpha12 * alpha23 * lam +
    2 * lam * alpha12 * alpha24 + 2  * alpha23 * alpha24 * lam + 2  * alpha24 * lam

mean(reaN)
mean(reaN^2) - mean(reaN)^2


## Simulation avec notre fonction
source("sample_MPMRF.R")

realisations2 <- rMPMRF(nn, adj, lam)

colMeans(realisations2)
apply(realisations2, 2, var)
cor(realisations2[, 1], realisations2[, 2])
cor(realisations2[, 2], realisations2[, 3])
cor(realisations2[, 2], realisations2[, 4])
cor(realisations2[, 1], realisations2[, 3])
cor(realisations2[, 1], realisations2[, 4])

reaN2 <- rowSums(realisations2)
mean(reaN2)
mean(reaN2^2) - mean(reaN2)^2


## pmf conjointe à la main et la nôtre
source("pdf_MPMRF.R")

fmp_n <- function(x1, x2, x3, x4, a12, a23, a24, lambda)
{
    dpois(x1, lambda) *
        (sum(dpois(x2 - 0:min(x1, x2), lambda * (1 - a12)) *
                 dbinom(0:min(x1, x2), x1, a12))) *
        (sum(dpois(x3 - 0:min(x2, x3), lambda * (1 - a23)) *
                 dbinom(0:min(x2, x3), x2, a23))) *
        (sum(dpois(x4 - 0:min(x2, x4), lambda * (1 - a24)) *
                 dbinom(0:min(x2, x4), x2, a24)))
}

xx <- c(0, 0, 0, 0)
fmp_n(xx[1], xx[2], xx[3], xx[4], alpha12, alpha23, alpha24, lam)
dpois(0, lam) * dpois(0, lam * (1 -alpha12)) * dpois(0, lam * (1 - alpha23)) * dpois(0, lam * (1 - alpha24))
pdf_MPMRF(adj, lam, xx, 1)

xx <- c(2, 5, 7, 8)
fmp_n(xx[1], xx[2], xx[3], xx[4], alpha12, alpha23, alpha24, lam)
pdf_MPMRF(adj, lam, xx, 1)


xx <- c(2, 5, 7, 8)
fmp_n(xx[1], xx[2], xx[3], xx[4], alpha12, alpha23, alpha24, lam)
pdf_MPMRF(adj, lam, xx, 1)


xx <- c(8, 8, 8, 5)
fmp_n(xx[1], xx[2], xx[3], xx[4], alpha12, alpha23, alpha24, lam)
pdf_MPMRF(adj, lam, xx, 1)


xx <- c(5, 5, 5, 5)
fmp_n(xx[1], xx[2], xx[3], xx[4], alpha12, alpha23, alpha24, lam)
pdf_MPMRF(adj, lam, xx, 1)

grid <- expand.grid(0:10, 0:10, 0:10, 0:10)
valid_pmf <- matrix(numeric(2 * nrow(grid)), ncol = 2)
for (i in 1:nrow(grid))
{
    (a <- fmp_n(grid[i, 1], grid[i, 2], grid[i, 3], grid[i, 4], alpha12, alpha23, alpha24, lam))
    (b <- pdf_MPMRF(adj, lam, unlist(grid[i, ]), 1))
    valid_pmf[i, ] <- c(a, b)
}
valid_pmf
colSums(valid_pmf)
decal <- apply(valid_pmf, 1, function(a) (abs(diff(a)) > 1e-18))

sum(decal)
valid_pmf[decal, ]

#### fgp M à la main

fgpN <- function(t1, t2, t3, t4, a12, a23, a24, lam)
{
    expo1 <- t1 * (1 - a12 + a12 * (t2 * (1 - a23 + a23 * t3) * (1 - a24 + a24 * t4))) - 1
    expo2 <- t2 * (1 - a23 + a23 * t3) * (1 - a24 + a24 * t4) - 1
    expo3 <- t3 - 1
    expo4 <- t4 - 1

    exp(lam * expo1) * exp(lam * (1 - a12) * expo2) *
        exp(lam * (1 - a23) * expo3) * exp(lam * (1 - a24) * expo4)
}

fgpN(1, 1, 1, 1, alpha12, alpha23, alpha24, lam)

nfft <- 2^15

fb <- c(0, 1, rep(0, nfft - 2))

ffb <- fft(fb)

ffM <- fgpN(ffb, ffb, ffb, ffb, alpha12, alpha23, alpha24, lam)
fM <- Re(fft(ffM, inverse = TRUE)) / nfft
sum(fM)

sum(0:(nfft - 1) * fM)
sum((0:(nfft - 1))^2 * fM) - sum(0:(nfft - 1) * fM)^2
4 * lam + 2 * alpha12 * lam + 2 * alpha23 * lam + 2 * alpha12 * alpha23 * lam +
    2 * lam * alpha12 * alpha24 + 2  * alpha23 * alpha24 * lam + 2  * alpha24 * lam


#### fgp notre fonction pour M
source("pdf_M.R")
fM_ <- pdf_M(adj, lam)

sum(fM_)
sum(0:(nfft - 1) * fM_)
sum((0:(nfft - 1))^2 * fM_) - sum(0:(nfft - 1) * fM_)^2
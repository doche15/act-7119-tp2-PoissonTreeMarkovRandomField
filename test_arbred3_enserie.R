###
### ACT-7119
### Test arbe en série de d = 3
### [Côté et al., 2025].
###
##

library(igraph)
source("graph_A.R")

alpha12 <- 0.4
alpha23 <- 0.5
A <- matrix(c(1, alpha12, 0,
              alpha12, 1, alpha23,
              0, alpha23, 1), ncol = 3, byrow = TRUE)

see_tree_graph(A, root_node = 1, titre = "Graphique en série d = 3")

lam <- 5

## Simulation à la main

rea <- function(a12, a23, lambda)
{
    N1 <- rpois(1, lambda)
    prop2 <- rbinom(1, N1, a12)
    inno2 <- rpois(1, lambda * (1 - a12))

    N2 <- prop2 + inno2

    prop3 <- rbinom(1, N2, a23)
    inno3 <- rpois(1, lambda * (1 - a23))

    N3 <- prop3 + inno3

    c(N1, N2, N3)
}

set.seed(20250301)
nn <- 100000
realisations <- t(replicate(nn, rea(alpha12, alpha23, lam)))
reaN <- rowSums(realisations)

colMeans(realisations)
apply(realisations, 2, var)
cor(realisations[, 1], realisations[, 2])
cor(realisations[, 2], realisations[, 3])

3 * lam + 2 * alpha12 * lam + 2 * alpha23 * lam + 2 * alpha12 * alpha23 * lam

mean(reaN)
mean(reaN^2) - mean(reaN)^2


## Simulation avec notre fonction
source("sample_MPMRF.R")

realisations2 <- rMPMRF(nn, A, lam)

colMeans(realisations2)
apply(realisations2, 2, var)
cor(realisations2[, 1], realisations2[, 2])
cor(realisations2[, 2], realisations2[, 3])

reaN2 <- rowSums(realisations2)
mean(reaN2)
mean(reaN2^2) - mean(reaN2)^2


## pmf conjointe à la main et la nôtre
source("pdf_MPMRF.R")

fmp_n <- function(x1, x2, x3, a12, a23, lambda)
{
    dpois(x1, lambda) *
        (sum(dpois(x2 - 0:min(x1, x2), lambda * (1 - a12)) *
                 dbinom(0:min(x1, x2), x1, a12))) *
        (sum(dpois(x3 - 0:min(x2, x3), lambda * (1 - a23)) *
                 dbinom(0:min(x2, x3), x2, a23)))
}

xx <- c(0, 0, 0)
fmp_n(xx[1], xx[2], xx[3], alpha12, alpha23, lam)
dpois(0, lam) * dpois(0, lam * (1 -alpha12)) * dpois(0, lam * (1 - alpha23))
pdf_MPMRF(A, lam, xx, 1)

xx <- c(2, 5, 7)
fmp_n(xx[1], xx[2], xx[3], alpha12, alpha23, lam)
pdf_MPMRF(A, lam, xx, 1)


xx <- c(8, 34, 5)
fmp_n(xx[1], xx[2], xx[3], alpha12, alpha23, lam)
pdf_MPMRF(A, lam, xx, 1)


xx <- c(8, 8, 8)
fmp_n(xx[1], xx[2], xx[3], alpha12, alpha23, lam)
pdf_MPMRF(A, lam, xx, 1)


xx <- c(5, 5, 5)
fmp_n(xx[1], xx[2], xx[3], alpha12, alpha23, lam)
pdf_MPMRF(A, lam, xx, 1)

grid <- expand.grid(0:20, 0:20, 0:20)
valid_pmf <- matrix(numeric(2 * nrow(grid)), ncol = 2)
for (i in 1:nrow(grid))
{
    (a <- fmp_n(grid[i, 1], grid[i, 2], grid[i, 3], alpha12, alpha23, lam))
    (b <- pdf_MPMRF(A, lam, unlist(grid[i, ]), 1))
    valid_pmf[i, ] <- c(a, b)
}
valid_pmf


#### fgp M à la main

fgpN <- function(t1, t2, t3, a12, a23, lam)
{
    expo1 <- t1 * (1 - a12 + a12 * (t2 * (1 - a23 + a23 * t3))) - 1
    expo2 <- t2 * (1 - a23 + a23 * t3) - 1
    expo3 <- t3 - 1

    exp(lam * expo1) * exp(lam * (1 - a12) * expo2) * exp(lam * (1 - a23) * expo3)
}

fgpN(1, 1, 1, alpha12, alpha23, lam)

nfft <- 2^15

fb <- c(0, 1, rep(0, nfft - 2))

ffb <- fft(fb)

ffM <- fgpN(ffb, ffb, ffb, alpha12, alpha23, lam)
fM <- Re(fft(ffM, inverse = TRUE)) / nfft
sum(fM)

sum(0:(nfft - 1) * fM)
sum((0:(nfft - 1))^2 * fM) - sum(0:(nfft - 1) * fM)^2
3 * lam + 2 * alpha12 * lam + 2 * alpha23 * lam + 2 * alpha12 * alpha23 * lam


#### fgp notre fonction pour M
source("pdf_M.R")
fM_ <- pdf_M(A, lam)

sum(fM_)
sum(0:(nfft - 1) * fM_)
sum((0:(nfft - 1))^2 * fM_) - sum(0:(nfft - 1) * fM_)^2

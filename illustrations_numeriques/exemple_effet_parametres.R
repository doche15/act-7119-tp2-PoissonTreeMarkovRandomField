###
### Travail 2, ACT-7119
### Exemple effet de la structure
###

library(igraph)

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

d <- 25
A <- diag(25)

lam <- 2

source("sample_MPMRF.R")
source("pdf_M.R")
source("graph_A.R")
source("algorithm_4.R")


#### Structure 3
## Arètes


## Cas 1 a)
alpha <- 0.3
beta <- 0.3
aretes <- matrix(c(1, 6, 6, 7, 6, 8, 6, 9, 1, 2, 2, 3, 2, 4, 2, 5, 1, 10, 10,
                   11, 10, 12, 10, 13, 1, 14, 14, 15, 14, 16, 14, 17, 1, 18, 18,
                   19, 18, 20, 18, 21, 1, 22, 22, 23, 22, 24, 22, 25),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, c(alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3)))
A3 <- A
for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A3[u, v] <- aretes_alpha[i, 3]
    A3[v, u] <- aretes_alpha[i, 3]
}

contrib_tvar <- function(fM, alloc, lambda, kappa)
{
    FM <- cumsum(fM)
    varr <- min(which(FM >= kappa)) - 1
    aa <- lambda - sum(alloc[0:varr + 1]) + (FM[varr + 1] - kappa) /fM[varr + 1] *
        alloc[varr + 1]

    1 / (1 - kappa) * aa
}

## Simulations
nn <- 5e6
set.seed(66)
realisations <- rMPMRF(nn, A3, lam)

reaM3 <- rowSums(realisations)

## Exact
nfft <- 2^15
fM3 <- pdf_M(A3, lam, nfft)

mm <- 0:(nfft - 1)

sum(mm * fM3)
mean(reaM3)


sum(mm^2 * fM3) - sum(mm * fM3)^2
mean(reaM3^2) - mean(reaM3)^2

cumsum(fM3)[20]
mean(reaM3 <= 19)


(v0.9 <- min(which(cumsum(fM3) >= 0.9)) - 1)
sort(reaM3)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM3) / (1 - 0.9) + v0.9
mean(sort(reaM3)[(0.9 * nn):nn])

alloc1 <- exp_alloc_Nv(A3, 1, lam, 1, nfft)
alloc6 <- exp_alloc_Nv(A3, 1, lam, 6, nfft)
alloc16 <- exp_alloc_Nv(A3, 1, lam, 16, nfft)

contrib_tvar(fM3, alloc1, lam, 0.9)
contrib_tvar(fM3, alloc6, lam, 0.9)
contrib_tvar(fM3, alloc16, lam, 0.9)


## Cas 1 b)
alpha <- 0.6
beta <- 0.6
aretes <- matrix(c(1, 6, 6, 7, 6, 8, 6, 9, 1, 2, 2, 3, 2, 4, 2, 5, 1, 10, 10,
                   11, 10, 12, 10, 13, 1, 14, 14, 15, 14, 16, 14, 17, 1, 18, 18,
                   19, 18, 20, 18, 21, 1, 22, 22, 23, 22, 24, 22, 25),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, c(alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3)))
A3 <- A
for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A3[u, v] <- aretes_alpha[i, 3]
    A3[v, u] <- aretes_alpha[i, 3]
}

contrib_tvar <- function(fM, alloc, lambda, kappa)
{
    FM <- cumsum(fM)
    varr <- min(which(FM >= kappa)) - 1
    aa <- lambda - sum(alloc[0:varr + 1]) + (FM[varr + 1] - kappa) /fM[varr + 1] *
        alloc[varr + 1]

    1 / (1 - kappa) * aa
}

## Simulations
nn <- 5e6
set.seed(66)
realisations <- rMPMRF(nn, A3, lam)

reaM3 <- rowSums(realisations)

## Exact
nfft <- 2^15
fM3 <- pdf_M(A3, lam, nfft)

mm <- 0:(nfft - 1)

sum(mm * fM3)
mean(reaM3)


sum(mm^2 * fM3) - sum(mm * fM3)^2
mean(reaM3^2) - mean(reaM3)^2

cumsum(fM3)[20]
mean(reaM3 <= 19)


(v0.9 <- min(which(cumsum(fM3) >= 0.9)) - 1)
sort(reaM3)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM3) / (1 - 0.9) + v0.9
mean(sort(reaM3)[(0.9 * nn):nn])

alloc1 <- exp_alloc_Nv(A3, 1, lam, 1, nfft)
alloc6 <- exp_alloc_Nv(A3, 1, lam, 6, nfft)
alloc16 <- exp_alloc_Nv(A3, 1, lam, 16, nfft)

contrib_tvar(fM3, alloc1, lam, 0.9)
contrib_tvar(fM3, alloc6, lam, 0.9)
contrib_tvar(fM3, alloc16, lam, 0.9)

## Cas 2 a)
alpha <- 0.5
beta <- 0.75
aretes <- matrix(c(1, 6, 6, 7, 6, 8, 6, 9, 1, 2, 2, 3, 2, 4, 2, 5, 1, 10, 10,
                   11, 10, 12, 10, 13, 1, 14, 14, 15, 14, 16, 14, 17, 1, 18, 18,
                   19, 18, 20, 18, 21, 1, 22, 22, 23, 22, 24, 22, 25),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, c(alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3)))
A3 <- A
for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A3[u, v] <- aretes_alpha[i, 3]
    A3[v, u] <- aretes_alpha[i, 3]
}

contrib_tvar <- function(fM, alloc, lambda, kappa)
{
    FM <- cumsum(fM)
    varr <- min(which(FM >= kappa)) - 1
    aa <- lambda - sum(alloc[0:varr + 1]) + (FM[varr + 1] - kappa) /fM[varr + 1] *
        alloc[varr + 1]

    1 / (1 - kappa) * aa
}

## Simulations
nn <- 5e6
set.seed(66)
realisations <- rMPMRF(nn, A3, lam)

reaM3 <- rowSums(realisations)

## Exact
nfft <- 2^15
fM3 <- pdf_M(A3, lam, nfft)

mm <- 0:(nfft - 1)

sum(mm * fM3)
mean(reaM3)


sum(mm^2 * fM3) - sum(mm * fM3)^2
mean(reaM3^2) - mean(reaM3)^2

cumsum(fM3)[20]
mean(reaM3 <= 19)


(v0.9 <- min(which(cumsum(fM3) >= 0.9)) - 1)
sort(reaM3)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM3) / (1 - 0.9) + v0.9
mean(sort(reaM3)[(0.9 * nn):nn])

alloc1 <- exp_alloc_Nv(A3, 1, lam, 1, nfft)
alloc6 <- exp_alloc_Nv(A3, 1, lam, 6, nfft)
alloc16 <- exp_alloc_Nv(A3, 1, lam, 16, nfft)

contrib_tvar(fM3, alloc1, lam, 0.9)
contrib_tvar(fM3, alloc6, lam, 0.9)
contrib_tvar(fM3, alloc16, lam, 0.9)

## Cas 2 b)
alpha <- 0.5
beta <- 0.9
aretes <- matrix(c(1, 6, 6, 7, 6, 8, 6, 9, 1, 2, 2, 3, 2, 4, 2, 5, 1, 10, 10,
                   11, 10, 12, 10, 13, 1, 14, 14, 15, 14, 16, 14, 17, 1, 18, 18,
                   19, 18, 20, 18, 21, 1, 22, 22, 23, 22, 24, 22, 25),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, c(alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3)))
A3 <- A
for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A3[u, v] <- aretes_alpha[i, 3]
    A3[v, u] <- aretes_alpha[i, 3]
}

contrib_tvar <- function(fM, alloc, lambda, kappa)
{
    FM <- cumsum(fM)
    varr <- min(which(FM >= kappa)) - 1
    aa <- lambda - sum(alloc[0:varr + 1]) + (FM[varr + 1] - kappa) /fM[varr + 1] *
        alloc[varr + 1]

    1 / (1 - kappa) * aa
}

## Simulations
nn <- 5e6
set.seed(66)
realisations <- rMPMRF(nn, A3, lam)

reaM3 <- rowSums(realisations)

## Exact
nfft <- 2^15
fM3 <- pdf_M(A3, lam, nfft)

mm <- 0:(nfft - 1)

sum(mm * fM3)
mean(reaM3)


sum(mm^2 * fM3) - sum(mm * fM3)^2
mean(reaM3^2) - mean(reaM3)^2

cumsum(fM3)[20]
mean(reaM3 <= 19)


(v0.9 <- min(which(cumsum(fM3) >= 0.9)) - 1)
sort(reaM3)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM3) / (1 - 0.9) + v0.9
mean(sort(reaM3)[(0.9 * nn):nn])

alloc1 <- exp_alloc_Nv(A3, 1, lam, 1, nfft)
alloc6 <- exp_alloc_Nv(A3, 1, lam, 6, nfft)
alloc16 <- exp_alloc_Nv(A3, 1, lam, 16, nfft)

contrib_tvar(fM3, alloc1, lam, 0.9)
contrib_tvar(fM3, alloc6, lam, 0.9)
contrib_tvar(fM3, alloc16, lam, 0.9)


## Cas 3a)
alpha <- 0.7
beta <- 0.4
aretes <- matrix(c(1, 6, 6, 7, 6, 8, 6, 9, 1, 2, 2, 3, 2, 4, 2, 5, 1, 10, 10,
                   11, 10, 12, 10, 13, 1, 14, 14, 15, 14, 16, 14, 17, 1, 18, 18,
                   19, 18, 20, 18, 21, 1, 22, 22, 23, 22, 24, 22, 25),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, c(alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3)))
A3 <- A
for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A3[u, v] <- aretes_alpha[i, 3]
    A3[v, u] <- aretes_alpha[i, 3]
}

contrib_tvar <- function(fM, alloc, lambda, kappa)
{
    FM <- cumsum(fM)
    varr <- min(which(FM >= kappa)) - 1
    aa <- lambda - sum(alloc[0:varr + 1]) + (FM[varr + 1] - kappa) /fM[varr + 1] *
        alloc[varr + 1]

    1 / (1 - kappa) * aa
}

## Simulations
nn <- 5e6
set.seed(66)
realisations <- rMPMRF(nn, A3, lam)

reaM3 <- rowSums(realisations)

## Exact
nfft <- 2^15
fM3 <- pdf_M(A3, lam, nfft)

mm <- 0:(nfft - 1)

sum(mm * fM3)
mean(reaM3)


sum(mm^2 * fM3) - sum(mm * fM3)^2
mean(reaM3^2) - mean(reaM3)^2

cumsum(fM3)[20]
mean(reaM3 <= 19)


(v0.9 <- min(which(cumsum(fM3) >= 0.9)) - 1)
sort(reaM3)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM3) / (1 - 0.9) + v0.9
mean(sort(reaM3)[(0.9 * nn):nn])

alloc1 <- exp_alloc_Nv(A3, 1, lam, 1, nfft)
alloc6 <- exp_alloc_Nv(A3, 1, lam, 6, nfft)
alloc16 <- exp_alloc_Nv(A3, 1, lam, 16, nfft)

contrib_tvar(fM3, alloc1, lam, 0.9)
contrib_tvar(fM3, alloc6, lam, 0.9)
contrib_tvar(fM3, alloc16, lam, 0.9)


## Cas 3b)
alpha <- 0.95
beta <- 0.4
aretes <- matrix(c(1, 6, 6, 7, 6, 8, 6, 9, 1, 2, 2, 3, 2, 4, 2, 5, 1, 10, 10,
                   11, 10, 12, 10, 13, 1, 14, 14, 15, 14, 16, 14, 17, 1, 18, 18,
                   19, 18, 20, 18, 21, 1, 22, 22, 23, 22, 24, 22, 25),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, c(alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3),
                                alpha, rep(beta, 3)))
A3 <- A
for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A3[u, v] <- aretes_alpha[i, 3]
    A3[v, u] <- aretes_alpha[i, 3]
}

contrib_tvar <- function(fM, alloc, lambda, kappa)
{
    FM <- cumsum(fM)
    varr <- min(which(FM >= kappa)) - 1
    aa <- lambda - sum(alloc[0:varr + 1]) + (FM[varr + 1] - kappa) /fM[varr + 1] *
        alloc[varr + 1]

    1 / (1 - kappa) * aa
}

## Simulations
nn <- 5e6
set.seed(66)
realisations <- rMPMRF(nn, A3, lam)

reaM3 <- rowSums(realisations)

## Exact
nfft <- 2^15
fM3 <- pdf_M(A3, lam, nfft)

mm <- 0:(nfft - 1)

sum(mm * fM3)
mean(reaM3)


sum(mm^2 * fM3) - sum(mm * fM3)^2
mean(reaM3^2) - mean(reaM3)^2

cumsum(fM3)[20]
mean(reaM3 <= 19)


(v0.9 <- min(which(cumsum(fM3) >= 0.9)) - 1)
sort(reaM3)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM3) / (1 - 0.9) + v0.9
mean(sort(reaM3)[(0.9 * nn):nn])

alloc1 <- exp_alloc_Nv(A3, 1, lam, 1, nfft)
alloc6 <- exp_alloc_Nv(A3, 1, lam, 6, nfft)
alloc16 <- exp_alloc_Nv(A3, 1, lam, 16, nfft)

contrib_tvar(fM3, alloc1, lam, 0.9)
contrib_tvar(fM3, alloc6, lam, 0.9)
contrib_tvar(fM3, alloc16, lam, 0.9)
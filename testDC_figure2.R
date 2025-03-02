library(igraph)
source("graph_A.R")
source("sample_MPMRF.R")
source("pdf_M.R")
source("pdf_MPMRF.R")

alpha12 <- 0.99
alpha13 <- 0.9
alpha34 <- 0.1453
alpha35 <- 0.1
alpha46 <- 0.1134
alpha47 <- 0.9
lam <- 6
adj <- matrix(c(1, alpha12, alpha13, 0, 0, 0, 0,
                alpha12, 1, 0, 0, 0, 0, 0,
                alpha13, 0, 1, alpha34, alpha35, 0, 0,
                0, 0, alpha34, 1, 0, alpha46, alpha47,
                0, 0, alpha35, 0, 1, 0, 0,
                0, 0, 0, alpha46, 0, 1, 0,
                0, 0, 0, alpha47, 0, 0, 1), ncol = 7, byrow = TRUE)

see_tree_graph(adj, 1, titre = "Figure 2")

nn <- 1000000
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
vv <- sum((0:(nfft - 1))^2 * fM) - sum((0:(nfft - 1)) * fM)^2
nn <- 50000000
realisations <- rMPMRF(nn, adj, lam)
reaN <- rowSums(realisations)
mean(reaN)
mean(reaN^2) - mean(reaN)^2; vv

cumsum(fM)[30 + 1]
mean(reaN <= 30)

cumsum(fM)[46 + 1]
mean(reaN <= 46)

# fgp de Jean-Gab
fgp_N_figure_2 <- function(t_, a12, a13, a34, a35, a46, a47, lam)
{
    eta7 <- t_[7]
    eta6 <- t_[6]
    eta5 <- t_[5]
    eta4 <- t_[4] * (1 - a46 + a46 * eta6) * (1 - a47 + a47 * eta7)
    eta3 <- t_[3] * (1 - a34 + a34 * eta4) * (1 - a35 + a35 * eta5)
    eta2 <- t_[2]
    eta1 <- t_[1] * (1 - a12 + a12 * eta2) * (1 - a13 + a13 * eta3)

    exp(lam * (1 - 0) * (eta1 - 1)) *
        exp(lam * (1 - a12) * (eta2 - 1)) *
        exp(lam * (1 - a13) * (eta3 - 1)) *
        exp(lam * (1 - a34) * (eta4 - 1)) *
        exp(lam * (1 - a35) * (eta5 - 1)) *
        exp(lam * (1 - a46) * (eta6 - 1)) *
        exp(lam * (1 - a47) * (eta7 - 1))
}

fgp_M_figure_2 <- function(t, a12, a13, a34, a35, a46, a47, lam)
{
    eta7 <- t
    eta6 <- t
    eta5 <- t
    eta4 <- t * (1 - a46 + a46 * eta6) * (1 - a47 + a47 * eta7)
    eta3 <- t * (1 - a34 + a34 * eta4) * (1 - a35 + a35 * eta5)
    eta2 <- t
    eta1 <- t * (1 - a12 + a12 * eta2) * (1 - a13 + a13 * eta3)

    # etas <- c(eta1 = eta1,
    #           eta2 = eta2,
    #           eta3 = eta3,
    #           eta4 = eta4,
    #           eta5 = eta5,
    #           eta6 = eta6,
    #           eta7 = eta7)
    # print(etas)

    exp(lam * (1 - 0) * (eta1 - 1)) *
        exp(lam * (1 - a12) * (eta2 - 1)) *
        exp(lam * (1 - a13) * (eta3 - 1)) *
        exp(lam * (1 - a34) * (eta4 - 1)) *
        exp(lam * (1 - a35) * (eta5 - 1)) *
        exp(lam * (1 - a46) * (eta6 - 1)) *
        exp(lam * (1 - a47) * (eta7 - 1))
}

ffb <- fft(c(0, 1, rep(0, nfft - 2)))

ffm <- fgp_M_figure_2(ffb, alpha12, alpha13, alpha34, alpha35, alpha46, alpha47, lam)
fM_maison <- Re(fft(ffm, inverse = TRUE)) / nfft
sum(fM_maison)

sum((0:(nfft - 1)) * fM_maison)


mean(reaN)
mean(reaN^2) - mean(reaN)^2; vv; sum((0:(nfft - 1))^2 * fM_maison) - sum((0:(nfft - 1)) * fM_maison)^2

cumsum(fM)[30 + 1]
cumsum(fM_maison)[30 + 1]
mean(reaN <= 30)

cumsum(fM)[46 + 1]
mean(reaN <= 46)


### Test avec la pmf conjointe

fmp_n <- function(x1, x2, x3, x4, x5, x6, x7, a12, a13, a34,
                  a35, a46, a47, lambda)
{
    dpois(x1, lambda) *
        (sum(dpois(x2 - 0:min(x1, x2), lambda * (1 - a12)) *
                 dbinom(0:min(x1, x2), x1, a12))) *
        (sum(dpois(x3 - 0:min(x1, x3), lambda * (1 - a13)) *
                 dbinom(0:min(x1, x3), x1, a13))) *
        (sum(dpois(x4 - 0:min(x3, x4), lambda * (1 - a34)) *
                 dbinom(0:min(x3, x4), x3, a34))) *
    (sum(dpois(x5 - 0:min(x3, x5), lambda * (1 - a35)) *
             dbinom(0:min(x3, x5), x3, a35))) *
        (sum(dpois(x6 - 0:min(x4, x6), lambda * (1 - a46)) *
                 dbinom(0:min(x4, x6), x4, a46))) *
        (sum(dpois(x7 - 0:min(x4, x7), lambda * (1 - a47)) *
                 dbinom(0:min(x4, x7), x4, a47)))
}

lam <- 6

fmp_n(0, 0, 0, 0, 0, 0, 0, alpha12, alpha13, alpha34, alpha35, alpha46, alpha47,
      lam)

dpois(0, lam) * dpois(0, lam * (1 - alpha12)) * dpois(0, lam * (1 - alpha13)) *
    dpois(0, lam * (1 - alpha35)) * dpois(0, lam * (1 - alpha34)) * dpois(0, lam * (1 - alpha46)) *
    dpois(0, lam * (1 - alpha47))

fmp_n(0, 0, 0, 0, 0, 0, 1, alpha12, alpha13, alpha34, alpha35, alpha46, alpha47,
      lam)

dpois(0, lam) * dpois(0, lam * (1 - alpha12)) * dpois(0, lam * (1 - alpha13)) *
    dpois(0, lam * (1 - alpha35)) * dpois(0, lam * (1 - alpha34)) * dpois(0, lam * (1 - alpha46)) *
    dpois(1, lam * (1 - alpha47))

pdf_MPMRF(adj, lam, c(rep(0, 6), 1), 1)

fmp_n(2, 3, 4, 0, 0, 1, 1, alpha12, alpha13, alpha34, alpha35, alpha46, alpha47,
      lam)
pdf_MPMRF(adj, lam, c(2, 3, 4, 0, 0, 1,  1), 1)


grid <- expand.grid(3:5, 3:5, 3:5, 1:5, 3:5, 3:5, 3:5)
valid_pmf <- matrix(numeric(2 * nrow(grid)), ncol = 2)
for (i in 1:nrow(grid))
{
    print(i)
    (a <- fmp_n(grid[i, 1], grid[i, 2], grid[i, 3], grid[i, 4],
                grid[i, 5], grid[i, 6], grid[i, 7], alpha12, alpha13, alpha34,
                alpha35, alpha46, alpha47, lam))
    (b <- pdf_MPMRF(adj, lam, unlist(grid[i, ]), 1))
    valid_pmf[i, ] <- c(a, b)
}
valid_pmf
colSums(valid_pmf)
decal <- apply(valid_pmf, 1, function(a) (abs(diff(a)) > 1e-18))
sum(decal) # :)



fmp_MPMRF <- function(x, A, lam)
{
    d <- nrow(A)

    prod_vec <- numeric(d - 1)

    for (i in 2:d)
    {
        pik <- min(which(A[i, ] > 0))

        prod_vec[i-1] <- sum(dpois(x[i] - 0:(min(x[pik], x[i])), lam *
                                       (1 - A[pik, i])) *
                dbinom(0:(min(x[pik], x[i])), x[pik], A[pik, i]))
    }
    dpois(x[1], lam) * prod(prod_vec)
}

fmp_n(2, 3, 4, 0, 0, 1, 1, alpha12, alpha13, alpha34, alpha35, alpha46, alpha47,
      lam)
pdf_MPMRF(adj, lam, c(2, 3, 4, 0, 0, 1,  1), 1)
fmp_MPMRF(c(2, 3, 4, 0, 0, 1,  1), adj, lam)

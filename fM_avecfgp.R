###
### Travail 2, ACT-7119
### Essai loi de M par la fgp
###

source("pgf_MPMRF.R")
source("pdf_MPMRF.R")
source("sample_MPMRF.R")
alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances
lam <- 1

# Matrice adjacente
A = matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

A = matrix(c(1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1),
           nrow = 4,
           byrow = TRUE)

d <- nrow(A)

## Va b dégénérée
nfft <- 2^9
fb <- numeric(nfft)
fb[2] <- 1

ffb <- fft(fb)


ffM <- sapply(1:nfft, function(a) pgf_MPMRF(A, lam, rep(ffb[a], d), root_node = 2))

ffM <- sapply(1:nfft, function(a) sum(rep(ffb[a], d)^a * pdf_MPMRF(A, lam, rep(a, d), 2)))

fM <- Re(fft(ffM, inverse = TRUE)) / nfft
sum(fM)

sum((0:(nfft - 1)) * fM)

## Variance théorique
d * lam + lam * (2 * alpha12 + 2 * alpha12 * 0 + 2 * alpha12 * alpha24 + 2 * alpha23 +
    2 * alpha24 + 2 * alpha23 * alpha24)

## Variance avec fM
sum((0:(nfft - 1))^2 * fM) - sum((0:(nfft - 1)) * fM)^2


rea <- rMPMRF(1000000, A, lam)

M <- rowSums(rea)

mean(M)
mean(M^2) - mean(M)^2


fM[1]
dpois(0, lam) * dpois(0, lam * (1-alpha12)) * dpois(0, lam * (1 - alpha24)) * dpois(0, lam * (1-alpha23))

pdf_MPMRF(A, lam, c(0, 0, 0, 0), 2)

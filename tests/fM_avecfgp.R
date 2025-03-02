###
### Travail 2, ACT-7119
### Essai loi de M par la fgp
###

source("pgf_MPMRF.R")
source("pdf_MPMRF.R")
source("sample_MPMRF.R")
source("pdf_M.R")
alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances
lam <- 1

# Matrice adjacente
A = matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

d <- nrow(A)

## Va b dégénérée
nfft <- 2^7
fb <- numeric(nfft)
fb[2] <- 1

ffb <- fft(fb)


ffM <- sapply(1:nfft, function(a) pgf_MPMRF(A, lam, rep(ffb[a], d), root_node = 1))

fM <- Re(fft(ffM, inverse = TRUE)) / nfft
sum(fM)

sum((0:(nfft - 1)) * fM)

## Variance théorique
d * lam + lam * (2 * alpha12 + 2 * alpha12 * 0 + 2 * alpha12 * alpha24 + 2 * alpha23 +
    2 * alpha24 + 2 * alpha23 * alpha24)

## Variance avec fM
sum((0:(20))^2 * fM[1:21]) - sum((0:(20)) * fM[1:21])^2


rea <- rMPMRF(1000000, A, lam)

M <- rowSums(rea)

mean(M)
mean(M^2) - mean(M)^2


fM[2]
dpois(0, lam) * dpois(0, lam * (1-alpha12)) * dpois(0, lam * (1 - alpha24)) * dpois(0, lam * (1-alpha23))

pdf_MPMRF(A, lam, c(1, 0, 0, 0), 1) + pdf_MPMRF(A, lam, c(0, 1, 0, 0), 1) + pdf_MPMRF(A, lam, c(0, 0, 1, 0), 1) + pdf_MPMRF(A, lam, c(0, 0, 0, 1), 1)


pdf_MPMRF(A, lam, c(1, 1, 0, 0), 1) + pdf_MPMRF(A, lam, c(0, 1, 1, 0), 1) + pdf_MPMRF(A, lam, c(0, 0, 1, 1), 1) +
    pdf_MPMRF(A, lam, c(0, 1, 0, 1), 1) + pdf_MPMRF(A, lam, c(1, 0, 1, 0), 1) + pdf_MPMRF(A, lam, c(1, 0, 0, 1), 1) +
    pdf_MPMRF(A, lam, c(2, 0, 0, 0), 1) + pdf_MPMRF(A, lam, c(0, 2, 0, 0), 1)  + pdf_MPMRF(A, lam, c(0, 0, 2, 0), 1)  +
    pdf_MPMRF(A, lam, c(0, 0, 0, 2), 1)



pdf_M(A, 1)[1:10]
fM[1:10]


grid <- expand.grid(0:50, 0:50, 0:50, 0:50)
M <- rowSums(grid)


M3 <- grid[M == 3, ]

sum(apply(M3, 1, function(a) pdf_MPMRF(A, lam, a, 1)))
fM[4]


m <- 20
Mm <- grid[M == m, ]
sum(apply(Mm, 1, function(a) pdf_MPMRF(A, lam, a, 1)))
fM[m+1]


### Test avec d = 2

# Matrice adjacente
A = matrix(c(1, alpha12,
             alpha12, 1),
           nrow = 2,
           byrow = TRUE)

pdf_MPMRF(A, lam, c(0, 1), 1)

dpois(0, lam) * dpois(1, (1 - alpha12) * lam) * dbinom(0, 0, alpha12)

nfft <- 2^7
fb <- numeric(nfft)
fb[2] <- 1

ffb <- fft(fb)
d <- nrow(A)

ffM <- sapply(1:nfft, function(a) pgf_MPMRF(A, lam, rep(ffb[a], d), root_node = 1))

fM <- Re(fft(ffM, inverse = TRUE)) / nfft
sum(fM)

sum((0:(nfft - 1)) * fM)
fM[1:10]

dpois(0, lam) * dpois(0, lam * (1 -alpha12))

pdf_MPMRF(A, lam, c(0, 1), root_node = 1) + pdf_MPMRF(A, lam, c(1, 0), root_node = 1)
fM[2]

grid <- expand.grid(0:20, 0:20)
M <- rowSums(grid)

for (i in 0:20)
{
    m <- i
    Mm <- grid[M == m, ]
    a <- sum(apply(Mm, 1, function(a) pdf_MPMRF(A, lam, a, 1))); b <- fM[m+1]
    print(c(a, b))
}

2 * lam + 2 * alpha12

sum((0:(nfft - 1))^2 * fM) - sum((0:(nfft - 1)) * fM)^2

rea <- rMPMRF(1000, A, lam)

M <- rowSums(rea)
var(M)



# Matrice adjacente
A = matrix(c(1, alpha12, 0,
             alpha12, 1, alpha23,
             0, alpha23, 1),
           nrow = 3,
           byrow = TRUE)

nfft <- 2^7
fb <- numeric(nfft)
fb[2] <- 1

ffb <- fft(fb)
d <- nrow(A)

ffM <- sapply(1:nfft, function(a) pgf_MPMRF(A, lam, rep(ffb[a], d), root_node = 1))

ffb[2]
pgf_MPMRF(A, lam, rep(ffb[2], d), root_node = 1)

fM <- Re(fft(ffM, inverse = TRUE)) / nfft
sum(fM)

sum((0:(nfft - 1)) * fM)
fM[1:10]

dpois(0, lam) * dpois(0, lam * (1 -alpha12)) * dpois(0, lam * (1 - alpha23))

pdf_MPMRF(A, lam, c(0, 1, 0), root_node = 1) + pdf_MPMRF(A, lam, c(1, 0, 0), root_node = 1) + + pdf_MPMRF(A, lam, c(0, 0, 1), root_node = 1)
fM[2]

grid <- expand.grid(0:20, 0:20, 0:20)
M <- rowSums(grid)

for (i in 0:20)
{
    m <- i
    Mm <- grid[M == m, ]
    a <- sum(apply(Mm, 1, function(a) pdf_MPMRF(A, lam, a, 1))); b <- fM[m+1]
    print(c(a, b))
}

3 * lam + 2 * alpha12 + 2 * alpha23

sum((0:(nfft - 1))^2 * fM) - sum((0:(nfft - 1)) * fM)^2

rea <- rMPMRF(1000, A, lam)

M <- rowSums(rea)
var(M)
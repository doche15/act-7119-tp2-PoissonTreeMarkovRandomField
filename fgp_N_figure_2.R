###
### Pmf de M
### Test en développant au long la fgp du vecteur aléatoire N de l'arbre de
### la Figure 2 de [Côté et al., 2025].
###
##

source("sample_MPMRF.R")

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


# Analyse
a <- c(0.1, 0.2, 0.3, 0.4, 0.3, 0.2)
lam <- 10

A <- matrix(c(1, a[1], a[2], 0, 0, 0, 0,
              a[1], 1, 0, 0, 0, 0, 0,
              a[2], 0, 1, a[3], a[4], 0, 0,
              0, 0, a[3], 1, 0, a[5], a[6],
              0, 0, a[4], 0, 1, 0, 0,
              0, 0, 0, a[5], 0, 1, 0,
              0, 0, 0, a[6], 0, 0, 1),
            nrow = 7)

nfft <- 2^18
n <- 1000000

# Échantillonnage
resultats_ech <- rMPMRF(n, A, lam)
resultats_ech_M <- rowSums(resultats_ech)
espM_emp <- mean(resultats_ech_M)
varM_emp <- var(resultats_ech_M)


# Méthode FFT avec fgp de N_
fmpB <- rep(0, nfft)
fmpB[2] <- 1

fftB <- fft(fmpB)

fftM <- fgp_M_figure_2(fftB, a[1], a[2], a[3], a[4], a[5], a[6], lam)

fM <- Re(fft(fftM, inverse = TRUE)) / nfft
fM <- pmax(fM, 0) # ajutement pour éviter les valeurs légèrement négatives
# qui pourraient être multipliées avec de grandes valeurs du support de M
sum(fM)

supportM <- 0:(nfft - 1)
espM_th <- sum(supportM * fM)
varM_th <- sum((supportM^2) * fM) - (sum(supportM * fM)^2)


c(espM_emp = espM_emp,
  espM_th = espM_th,
  varM_emp = varM_emp,
  varM_th = varM_th)


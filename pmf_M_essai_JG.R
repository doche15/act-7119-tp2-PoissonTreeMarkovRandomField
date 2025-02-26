###
### Travail 2, ACT-7119
### Algorithme 3 de [Côté et al., 2025]
###
##

calculate_pmf_M <- function(A, lambda, nfft = 2^17)
{
    # Vecteur b
    b <- rep(0, nfft)
    b[2] <- 1

    fft_b <- fft(b)
    fft_M <- numeric(nfft)

    dimension <- nrow(A)

    for (l in seq(nfft))
    {
        H <- matrix(rep(1, dimension^2), nrow = dimension)

        h_k <- numeric(dimension)

        pi_k <- numeric(dimension)

        for (k in rev(seq(dimension))[-dimension])
        {
            pi_k[abs(k - dimension) + 2] <- which(A[k,] > 0)[1]

            h_k[abs(k - dimension) + 2] <- fft_b[l] * prod(H[k,])

            H[pi_k[abs(k - dimension) + 2], k] <- (1 - A[pi_k[abs(k - dimension) + 2], k]) +
                A[pi_k[abs(k - dimension) + 2], k] *
                h_k[abs(k - dimension) + 2]
        }

        h_k[1] <- fft_b[l] * prod(H[1,])

        # print(prod(exp(lambda * (1 - A[pi_k,]) * (h_k - 1))))
        # print(unlist(sapply(seq(dimension), function(k) exp(lambda * (1 - A[pi_k[k], k]) * (h_k[k] - 1)))))

        fft_M[l] <- prod(unlist(sapply(seq(dimension), function(k) exp(lambda * (1 - A[pi_k[k], k]) * (h_k[k] - 1)))))
        # fft_M[l] <- prod(exp(lambda * (1 - A[pi_k,]) * (h_k - 1)))
    }

    Re(fft(fft_M, inverse = TRUE)) / nfft
}

# Même matrice A que la matrice A des autres codes.
A <- matrix(c(1, alpha12, 0, 0,
              alpha12, 1, alpha23, alpha24,
              0, alpha23, 1, 0,
              0, alpha24, 0, 1),
            nrow = 4,
            byrow = TRUE)

xx <- calculate_pmf_M(A, 1)
sum(xx)

xx[1] # devrait être égal à la valeur suivante.
dpois(0, 1 * (1 - 0)) *
    dpois(0, 1 * (1 - alpha12)) *
    dpois(0, 1 * (1 - alpha23)) *
    dpois(0, 1 * (1 - alpha24))

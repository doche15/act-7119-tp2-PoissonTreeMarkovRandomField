###
### Travail 2, ACT-7119
### Theorem 7 : Distribution of the Sum
###

pdf_M = function(A, lambda, nfft){
    # A : matrice d'adjacence
    # lambda : paramètre des lois de Poisson

    d = nrow(A)

    nfft = nfft

    b = c(0, 1, rep(0, nfft - 2)) # fmp v.a. dégénérée à 1

    phib = fft(b)

    phiM = numeric(nfft)

    stock_pi_k = numeric(d - 1)
    stock_h_k = numeric(d - 1)

    k_domain = rev(seq(d)[-1])

    for (l in seq(nfft)){

        H = matrix(data = 1,
                   nrow = d,
                   ncol = d) # all-1 matrix

        for (k in k_domain){

            pi_k = min(which(A[k,] > 0)) # infimum (devrait jamais donner infini parce
            # que sinon ça veut dire qu'il y a un node indépendant de tous les autres)
            # print(pi_k)

            h_k = phib[l] * prod(H[k,])

            H[pi_k, k] = 1 - A[pi_k, k] + A[pi_k, k] * h_k # overwrite H

            stock_pi_k[k - 1] = pi_k
            stock_h_k[k - 1] = h_k

        }

        h_1 = phib[l] * prod(H[1,])

        prod_vec = numeric(length(k_domain))

        for (k in k_domain){

            prod_vec[k - 1] = exp(lambda * (1 - A[stock_pi_k[k - 1], k]) *
                                      (stock_h_k[k - 1] - 1))

        }

        phiM[l] = prod(prod_vec) *
            exp(lambda * (h_1 - 1)) # ajouter k = 1

    }

    # output pdf_M
    Re(fft(phiM, inverse = T))/nfft

}

nfft = 2^16
# Validation pdf_M ----
alpha12 = 0.4 ; alpha23 = 0.4; alpha24 = 0.4 # dépendances
alpha34 <- 0.5

# Matrice adjacente
A = matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

A <- matrix(c(1, alpha12, 0,
              alpha12, 1, alpha23,
              0, alpha23, 1), nrow = 3)

A <- matrix(c(1, alpha12, 0, 0,
              alpha12, 1, alpha23, 0,
              0, alpha23, 1, alpha34,
              0, 0, alpha34, 1), nrow = 4)



a12 <- 0
a13 <- 0
a34 <- 0
a35 <- 0
a46 <- 0
a47 <- 0
A <- matrix(c(1, a12, a13, 0, 0, 0, 0,
              a12, 1, 0, 0, 0, 0, 0,
              a13, 0, 1, a34, a35, 0, 0,
              0, 0, a34, 1, 0, a46, a47,
              0, 0, a35, 0, 1, 0, 0,
              0, 0, 0, a46, 0, 1, 0,
              0, 0, 0, a47, 0, 0, 1),
            nrow = 7)
min(which(A[2,] > 0))
min(which(A[3,] > 0))
min(which(A[4,] > 0))




lambda = 10 # paramètre des lois de Poisson

fM = pdf_M(A, lambda, nfft)
sum(fM) # somme à 1

m = seq(nfft) - 1
esp_M = sum(m*fM) # changer les dépendances ne devrait pas avoir d'impact
var_M = sum((m^2)*fM) - esp_M^2 # augmenter la dépendance devrait augmenter la variance
print(esp_M)
print(var_M)

fM[1]
dpois(0, lambda * (1 - 0)) *
    dpois(0, lambda * (1 - a12)) *
    dpois(0, lambda * (1 - a13)) *
    dpois(0, lambda * (1 - a34)) *
    dpois(0, lambda * (1 - a35)) *
    dpois(0, lambda * (1 - a46)) *
    dpois(0, lambda * (1 - a47))

xx <- rMPMRF(10000000, A, lambda)
Mm <- rowSums(xx)

c(espth = esp_M, espemp = mean(Mm))
c(esp2th = sum((m^2)*fM), esp2emp = mean(Mm^2))
c(varth = var_M, varemp = var(Mm))

mean((sapply(0:60, function(x) mean(Mm == x)) - fM[1:61]) / fM[1:61])
(sapply(0:60, function(x) mean(Mm == x)) - fM[1:61]) / fM[1:61]

comparaison <- matrix(c(fM[1:61],
                        sapply(0:60, function(x) mean(Mm == x)),
                        100 * (sapply(0:60, function(x) mean(Mm == x)) - fM[1:61]) / fM[1:61]), ncol = 3, byrow = FALSE)
colnames(comparaison) <- c("TH", "EMP", "ECART (%)")

comparaison
# mean(Mm)
# var(Mm)
#
# mean(Mm == 0)

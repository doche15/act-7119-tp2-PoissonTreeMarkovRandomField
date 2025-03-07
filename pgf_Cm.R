###
### Travail 2, ACT-7119
### Theorem 7 : Distribution of the Sum
###
###     PGF of the secondary distribution
###

library(igraph)
source("find_parent_of_v.R")
source("nu_v_Tr.R")
source("pdf_MPMRF.R")
source("pdf_M.R")
pgf_Cm <-  function(A, lambda, t, root_node){
    # A : matrice d'adjacence
    # lambda : paramètre des lois de Poisson
    # t
    # root_node : numéro de la racine

    d <-  nrow(A)

    sum_vec <-  numeric(d)

    alpha_sum <- sum(A[lower.tri(A)])

    i <-  0
    for (v in seq(d)){

        i <-  i + 1

        alpha <-  ifelse(v == root_node, 0, A[find_parent_of_v(A, v, root_node), v])

        sum_vec[i] <-  (1 - alpha) / (d - alpha_sum) * nu(A, v, rep(t, d), root_node)

    }

    # output pgf
    #print("ok")
    sum(sum_vec)

}
# #
# # # Validation pgf_MPMRF ----
# alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances
# lam <- 3
# # # Matrice adjacente
# A = matrix(c(1, alpha12, 0, 0,
#              alpha12, 1, alpha23, alpha24,
#              0, alpha23, 1, 0,
#              0, alpha24, 0, 1),
#            nrow = 4,
#            byrow = TRUE)
#
# pgf_Cm(A, lam, 1, 1)
#
# alphas <- sum(A[lower.tri(A)])
# lamM <- lam * (nrow(A) - alphas)
#
# nfft <- 2^8
# fb <- c(0, 1, rep(0, nfft-2))
# ffb <- fft(fb)
#
# fM_aveccm <- Re(fft(exp(lamM * (sapply(ffb, function(a) pgf_Cm(A, lam, a, 1)) - 1)),
#                     inverse = TRUE))/nfft
#
# sum(fM_aveccm)
#
# fM <- pdf_M(A, lam, nfft = 2^8)
#
# cbind(fM_aveccm, fM) # ça marche!
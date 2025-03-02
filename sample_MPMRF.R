###
### Travail 2, ACT-7119
### Algorithm 2: Stochastic representation sampling method [Côté et al. 2025]
###

rMPMRF <-  function(n, A, lambda)
{
  # A : matrice d'adjacence

  d <-  nrow(A)

  # Matrice de réalisations.
  # n réalisations (lignes)
  # vecteur de longueur d (colonnes)
  N <- matrix(numeric(n * d), ncol = d)

  N[, 1] <- rpois(n, lambda)

  # Remplissage de la matrice par colonne
  for (k in 2:d)
  {
      pik <- min(which(A[k, ] > 0))
      Bk <- rbinom(n, N[, pik], A[pik, k])
      Lk <- rpois(n, lambda * (1 - A[pik, k]))
      N[, k] <- Bk + Lk
  }

  N
}
#
# # Validation pdf_MPMRF ----
# alpha12 <-  0.2 ; alpha23 <-  0.4; alpha24 <-  0.7 # dépendances
#
# # Matrice adjacente
# A <-  matrix(c(1, alpha12, 0, 0,
#              alpha12, 1, alpha23, alpha24,
#              0, alpha23, 1, 0,
#              0, alpha24, 0, 1),
#            nrow = 4,
#            byrow = TRUE)
#
# rea <- rMPMRF(100000, A, lambda = 4)
#
# colMeans(rea)
#
# cor(rea[, 1], rea[, 2])
# cor(rea[, 3], rea[, 2])
# cor(rea[, 3], rea[, 1])

###
### Travail 2, ACT-7119
### Algorithm 1 : Computing the variance-covariance matrix of N.
###
library(igraph)
# Calcul de A^(d) ----
compute_Ad = function(A){
  # A : matrice adjacente

  d = nrow(A)

  Aw = A

  for (w in seq(d)[-1]){

    A_temp = Aw

  # max-product matrix product
  for (i in seq(d)){

    for (j in seq(d)){

      vec = numeric(d)

      for (k in seq(d)){

        vec[k] = A_temp[i, k] * A[k, j]

      }

      Aw[i, j] = max(vec)

    }
   }
  }

 Aw

}

# Validation compute_Ad ----
alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances

# Matrice adjacente
A = matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

# A^(d)
compute_Ad(A)
# rho(N1, N3) = alpha12 * alpha23
alpha12 * alpha23
compute_Ad(A)[1, 3] # même résultat

# Calcul lambda*A^(d) ----
compute_var_cov_N = function(A, lambda){
  # A : matrice adjacente
  # lambda : paramètre des lois de Poisson

  lambda * compute_Ad(A)
}

# Validation compute_var_cov_N ----
# lambda*A^(d)
lambda = 3 # paramètre des lois de Poisson
compute_var_cov_N(A, lambda)

# cov(N1, N3) = lambda * alpha12 * alpha23
lambda * alpha12 * alpha23
compute_var_cov_N(A, lambda)[1, 3] # même résultat






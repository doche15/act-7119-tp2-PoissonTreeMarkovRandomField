###
### Travail 2, ACT-7119
### Theorem 3 : Joint Probability Mass Function.
###
library(igraph)
source("find_parent_of_v.R")

pdf_MPMRF = function(A, lambda, x_vec, root_node){
  # A : matrice d'adjacence
  # lambda : paramètre des lois de Poisson
  # x_vec : (x_1, ..., x_d)
  # root_node : numéro de la racine

  d = nrow(A)

  v_no_root = seq(d)[-root_node]

  prod_vec = numeric(d - 1)

  i = 0
  #browser()
  for (v in v_no_root){

    i = i + 1 # stock prod_vec

    # parent de v
    pav = find_parent_of_v(A, v, root_node)

    sum_v = 0

    for (k in seq(0, min(x_vec[pav], x_vec[v]))){

      #sum_v = sum_v +
      #(exp(-lambda * (1 - A[pav, v])) *
      #  (lambda * (1 - A[pav, v]))^(x_vec[v] - k)) /
      #  factorial(x_vec[v] - k) *
      #  choose(x_vec[pav], k) *
      #  A[pav, v]^k *
      #  (1 - A[pav, v])^(x_vec[pav] - k)

      ## Proposition pour simplifier :
      sum_v = sum_v + dpois(x_vec[v] - k, lambda *  (1 - A[pav, v])) *
                        dbinom(k, x_vec[pav], A[pav, v])

    }

    prod_vec[i] = sum_v
  }

  # output pdf
  dpois(x_vec[root_node], lambda) * prod(prod_vec)

}

### Version sans root_note
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

# # Validation pdf_MPMRF ----
# alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances
#
# # Matrice adjacente
# A = matrix(c(1, alpha12, 0, 0,
#              alpha12, 1, alpha23, alpha24,
#              0, alpha23, 1, 0,
#              0, alpha24, 0, 1),
#            nrow = 4,
#            byrow = TRUE)

# xi = 0:10
# grid = expand.grid(xi, xi, xi, xi)
# sum(apply(grid, 1, function(x) pdf_MPMRF(A, 1, x, 1))) # somme à 1
###
### Travail 2, ACT-7119
### Theorem 4 : Joint Probability Generating Function.
###
library(igraph)
source("find_parent_of_v.R")
source("nu_v_Tr.R")
source("pdf_MPMRF.R")

pgf_MPMRF = function(A, lambda, t_vec, root_node){
  # A : matrice adjacente
  # lambda : paramètre des lois de Poisson
  # t_vec : (t_1, ..., t_d)
  # root_node : numéro de la racine

  d = nrow(A)

  prod_vec = numeric(d)

  i = 0
  for (v in seq(d)){

    i = i + 1

    alpha = ifelse(v == root_node, 0, A[find_parent_of_v(A, v, root_node), v])

    prod_vec[i] = exp(lambda *
                        (1 - alpha) *
                        (nu(A, v, t_vec, root_node) - 1))

  }

  # output pgf
  prod(prod_vec)

}

# Validation pgf_MPMRF ----
alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances

# Matrice adjacente
A = matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

ti = rep(0, 4)
pgf_MPMRF(A, 1, ti, 1) # densité à 0
pdf_MPMRF(A, 1, ti, 1) # même chose

ti = rep(1, 4)
pgf_MPMRF(A, 1, ti, 1) # doit donner 1


ti <- c(0.3, 0.5, 0.6, 0.2)
pgf_MPMRF(A, 1, ti, 1)
pgf_MPMRF(A, 1, ti, 2)
pgf_MPMRF(A, 1, ti, 3)
pgf_MPMRF(A, 1, ti, 4)

###
### Travail 2, ACT-7119
### Algorithm 5 : Changing root in accordance to with an adjacency matrix is in topological order.
###

library(igraph)
source("graph_A.R")

reroot = function(A, root_node, new_root){
  # A : matrice d'adjacence
  # root_node : racine initiale
  # new_root : nouvelle racine voulue

  d = nrow(A)

  omega = new_root

  omega_vec = numeric()

  while (omega != root_node){

    omega_vec = c(omega_vec, omega)

    omega = min(which(A[omega,] > 0))

  }

  for (i in omega_vec){

    for (j in rev(seq(i - 1))){

      row_j = A[j, ] # ligne j
      row_jj = A[j + 1, ] # ligne j + 1

      # Switch
      A[j, ] = row_jj
      A[j + 1, ] = row_j

      col_j = A[, j] # colonne j
      col_jj = A[, j + 1] # colonne j + 1

      # Switch
      A[, j] = col_jj
      A[, j + 1] = col_j

    }

  }

  A

}

#
# # Validation graphique ----
# alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances
#
# # Matrice d'adjacence initiale
# A = matrix(c(1, alpha12, 0, 0,
#              alpha12, 1, alpha23, alpha24,
#              0, alpha23, 1, 0,
#              0, alpha24, 0, 1),
#            nrow = 4,
#            byrow = TRUE)
#
# see_tree_graph(A,
#                root_node = 1)
#
# # Matrice d'adjacence avec nouveau choix de racine
# A_prime = reroot(A,
#                  root_node = 1,
#                  new_root = 2)
#
# see_tree_graph(A_prime,
#                root_node = 1)
#
# Matrice d'adjacence avec nouveau choix de racine (problème)
# A_prime_prime = reroot(A,
#                        root_node = 1,
#                        new_root = 3)
# 
# see_tree_graph(A_prime_prime,
#                root_node = 1)


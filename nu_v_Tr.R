###
### Travail 2, ACT-7119
### Fonction récursive de nu_v^Tr(t_vdsc(v))
###
library(igraph)
source("find_ch_of_v.R")

# pas testée encore
nu = function(A, v, t_vec, root_node){
  # A : matrice adjacente
  # v : numéro du noeud cible
  # t_vec : (t_1, ..., t_d)
  # root_node : numéro de la racine

  ch_v = find_ch_of_v(A, v, root_node)

  if (length(ch_v) != 0){

  prod_vec = numeric(length(ch_v))

  i = 0
  for (j in ch_v){

    i = i + 1
    prod_vec[i] = 1 - A[v, j] + A[v, j] * nu(A, j, t_vec, root_node)

    }

  t_vec[v] * prod(prod_vec)

  } else {

    t_vec[v]

  }

}

# Validation nu ----
alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances

# Matrice adjacente (star)
A = matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

see_tree_graph(A, 
               root_node = 1,
               titre = "Arbre star")

ti = c(4, 7, 6, 2)

nu(A,
   v = 3,
   t_vec = ti,
   root_node = 1) # devrait être 6

ti[2] * (1 - A[2, 3] + A[2, 3] * (ti[3])) * (1 - A[2, 4] + A[2, 4] * (ti[4])) # valeur théorique
nu(A,
   v = 2,
   t_vec = ti,
   root_node = 1) # même chose

# Matrice adjacente (ligne)
alpha12 = 0.2 ; alpha23 = 0.4 # dépendances

A = matrix(c(1, alpha12, 0,
             alpha12, 1, alpha23,
             0, alpha23, 1),
           nrow = 3)

see_tree_graph(A, 
               root_node = 1,
               titre = "Arbre ligne")

ti = c(4, 7, 6)

nu(A,
   v = 3,
   t_vec = ti,
   root_node = 1) # devrait être 6

ti[2] * (1 - A[2, 3] + A[2, 3] * (ti[3])) # valeur théorique
nu(A,
   v = 2,
   t_vec = ti,
   root_node = 1)# même chose

# Matrice adjacente (hybride)
alpha12 = 0.2 ; alpha13 = 0.4 ; alpha34 = 0.7; alpha35 = 0.5# dépendances

A = matrix(c(1, alpha12, alpha13, 0, 0,
             alpha12, 1, 0, 0, 0,
             alpha13, 0, 1, alpha34, alpha35,
             0, 0, alpha34, 1, 0,
             0, 0, alpha35, 0, 1),
           nrow = 5)

see_tree_graph(A,
               root_node = 1,
               titre = "Arbre hybride")

ti = c(4, 7, 6, 5, 2)

nu(A,
   v = 4,
   t_vec = ti,
   root_node = 1) # devrait être 5

ti[1] * (1 - A[1, 2] + A[1, 2] * ti[2]) * (1 - A[1, 3] + A[1, 3] * ti[3] * (1 - A[3, 4] + A[3, 4] * ti[4]) * (1 - A[3, 5] + A[3, 5] * ti[5])) # valeur théorique
nu(A,
   v = 1,
   t_vec = ti,
   root_node = 1)# même chose


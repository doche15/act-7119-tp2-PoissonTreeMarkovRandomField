###
### Travail 2, ACT-7119
### Fonction pour identifier les enfants selon la structure
###
library(igraph)
source("find_parent_of_v.R")
source("graph_A.R")

# fonctionne pour le moment; à tester davantage
find_ch_of_v = function(A, v, root_node){
  # A : matrice adjacente
  # v : numéro de celui dont on veut les descendants
  # root_node : numéro de la racine

  d = nrow(A)

  ch_v = numeric()

  for (node in seq(d)[-c(v, root_node)]){

    pa_node = find_parent_of_v(A,
                               v = node,
                               root_node = root_node)

  if (pa_node == v){

    ch_v = c(ch_v, node)

    }

  }

  # output ch(v)
  ch_v

}
#
# # Validation find_ch_of_v ----
# alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances
#
# # Matrice adjacente (star)
# A = matrix(c(1, alpha12, 0, 0,
#              alpha12, 1, alpha23, alpha24,
#              0, alpha23, 1, 0,
#              0, alpha24, 0, 1),
#            nrow = 4,
#            byrow = TRUE)
#
# # ch(1) avec r = 1
# find_ch_of_v(A,
#              v = 1,
#              root_node = 1) # devrait être 2
#
# # ch(2) avec r = 1
# find_ch_of_v(A,
#              v = 2,
#              root_node = 1) # devrait être 3, 4
#
# # ch(3) avec r = 1
# find_ch_of_v(A,
#              v = 3,
#              root_node = 1) # aucun
#
# # ch(2) avec r = 2
# find_ch_of_v(A,
#              v = 2,
#              root_node = 2) # devrait être 1, 3, 4
#
# # Matrice adjacente (ligne)
# alpha12 = 0.2 ; alpha23 = 0.4 # dépendances
#
# A = matrix(c(1, alpha12, 0,
#              alpha12, 1, alpha23,
#              0, alpha23, 1),
#            nrow = 3)
#
# # ch(1) avec r = 1
# find_ch_of_v(A,
#              v = 1,
#              root_node = 1) # devrait être 2
#
# # ch(2) avec r = 1
# find_ch_of_v(A,
#              v = 2,
#              root_node = 1) # devrait être 3, 4
#
# # ch(3) avec r = 1
# find_ch_of_v(A,
#              v = 3,
#              root_node = 1) # aucun
#
# # ch(2) avec r = 2
# find_ch_of_v(A,
#              v = 2,
#              root_node = 2) # devrait être 1, 3, 4
#
# # Matrice adjacente (hybride)
# alpha12 = 0.2 ; alpha13 = 0.4 ; alpha34 = 0.7; alpha35 = 0.5# dépendances
#
# A = matrix(c(1, alpha12, alpha13, 0, 0,
#              alpha12, 1, 0, 0, 0,
#              alpha13, 0, 1, alpha34, alpha35,
#              0, 0, alpha34, 1, 0,
#              0, 0, alpha35, 0, 1),
#            nrow = 5)
#
# see_tree_graph(A,
#                root_node = 1,
#                titre = "Arbre hybride")
#
# # ch(2) avec r = 1
# find_ch_of_v(A,
#              v = 2,
#              root_node = 1) # aucun
#
# # ch(3) avec r = 1
# find_ch_of_v(A,
#              v = 3,
#              root_node = 1) # devrait être 4 et 5
#
# # ch(1) avec r = 1
# find_ch_of_v(A,
#              v = 1,
#              root_node = 1) # devrait être 2 et 3
#
# # change root pour r = 3
# see_tree_graph(A,
#                root_node = 3,
#                titre = "Arbre hybride")
#
# # ch(1) avec r = 3
# find_ch_of_v(A,
#              v = 1,
#              root_node = 3) # devrait être 2
#
# # ch(4) avec r = 3
# find_ch_of_v(A,
#              v = 4,
#              root_node = 3) # aucun
#
# # ch(3) avec r = 3
# find_ch_of_v(A,
#              v = 3,
#              root_node = 3) # devrait être 1, 4 et 5
#
#

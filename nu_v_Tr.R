###
### Travail 2, ACT-7119
### Fonction récursive de nu_v^Tr(t_vdsc(v))
###

source("find_dsc_of_v.R")

# pas testée encore
nu = function(A, v, t_vec, rooot_node){
  # A : matrice adjacente
  # v : numéro du noeud cible
  # t_vec : (t_1, ..., t_d)
  # root_node : numéro de la racine
  
  dsc_v = find_dsc_of_v(v)
  
  prod_vec = numeric(length(dsc_v))
  
  i = 0
  for (j in dsc_v){
    
    i = i + 1
    prod_vec[i] = 1 - A[v, j] + A[v, j] * nu(A, j, t_vec, root_node)
    
  }
  
  t_vec[v] * prod(prod_vec)
  
}





###
### Travail 2, ACT-7119
### Fonction récursive de nu_v^Tr(t_vdsc(v))
###

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

# Matrice adjacente
A = matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

ti = c(4,7,6,2)
nu(A, 3, ti, 1) # devrait être 6

7 * (1 - A[2, 3] + A[2, 3] * (6)) * (1 - A[2, 4] + A[2, 4] * (2)) # valeur théorique
nu(A, 2, ti, 1) # même chose




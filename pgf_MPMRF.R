###
### Travail 2, ACT-7119
### Theorem 4 : Joint Probability Generating Function.
###

source("find_parent_of_v.R")
source("nu_v_Tr.R")

# pas testée encore
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
    
    prod_vec[i] = exp(lambda * 
                        (1 - A[find_parent_of_v(v), v]) *
                        (nu(A, v, t_vec, rooot_node) - 1))
  
  }
  
  # output pgf
  prod(prod_vec)

}



###
### Travail 2, ACT-7119
### Fonction pour identifier les descendants selon la structure
###

source("find_parent_of_v.R")

# fonctionne pour le moment; à tester davantage
find_dsc_of_v = function(A, v, root_node){
  # A : matrice adjacente
  # v : numéro de celui dont on veut les descendants
  # root_node : numéro de la racine
  
  d = nrow(A)
  
  if (v != root_node){
  
  dsc_v = numeric()
  
  for (node in seq(d)[-c(v, root_node)]){
    
    pa_node = find_parent_of_v(A,
                               v = node,
                               root_node = root_node)
    
    if (pa_node == v){
      
      dsc_v = c(dsc_v, node)
      
    }
    
  }
  
  if (length(dsc_v) > 0){
    
    dsc_v = c(dsc_v, sapply(dsc_v, function(v) find_dsc_of_v(A, v, root_node)))
    
  }
  
  # output dsc(v)
  unlist(dsc_v)
  
  } else {
    
    # output dsc(v)
    seq(d)[-root_node]
    
  } 
  
}

# Validation find_dsc_of_v ----
alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances

# Matrice adjacente
A = matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

# dsc(1) avec r = 1
find_dsc_of_v(A,
              v = 1,
              root_node = 1) # devrait être 2, 3, 4

# dsc(2) avec r = 1
find_dsc_of_v(A,
              v = 2,
              root_node = 1) # devrait être 3, 4

# dsc(3) avec r = 1
find_dsc_of_v(A,
              v = 3,
              root_node = 1) # aucun

# dsc(2) avec r = 2
find_dsc_of_v(A,
              v = 2,
              root_node = 2) # devrait être 1, 3, 4

###
### Travail 2, ACT-7119
### Algorithm 5 : Computing expected allocations of Nv to M.
###

source("algorithm_5.R")

# à tester (devrait fonctionner quand pdf_M.R sera fonctionnel)
exp_alloc_Nv = function(A, root_node, lambda, v){
  # A : matrice d'adjacence
  # root_node : numéro de la racine initiale de A
  # lambda : paramètre des lois de Poisson
  # v : numéro du noeud dont on veut l'allocation
  
  A_prime = reroot(A, 
                   root_node = root_node, 
                   new_root = v)
  
  # Algorithm 3
  d = nrow(A)
  
  nfft = 2^15
  
  b = c(0, 1, rep(0, nfft - 2))
  
  phib = fft(b)
  
  phiM = numeric(nfft)
  
  stock_pi_k = numeric(d - 1)
  stock_h_k = numeric(d - 1)
  
  k_domain = rev(seq(d)[-1])
  
  for (l in seq(nfft)){
    
    H = matrix(data = 1,
               nrow = d,
               ncol = d)
    
    for (k in k_domain){
      
      pi_k = min(which(A[k,] > 0))
      
      h_k = phib[l] * prod(H[k,])
      
      H[pi_k, k] = 1 - A[pi_k, k] + A[pi_k, k] * h_k
      
      stock_pi_k[k - 1] = pi_k
      stock_h_k[k - 1] = h_k
      
    }
    
    h_1 = phib[l] * prod(H[1,])
    
    prod_vec = numeric(length(k_domain))
    
    for (k in k_domain){
      
      prod_vec[k - 1] = exp(lambda * (1 - A[stock_pi_k[k - 1], k]) * 
                              (stock_h_k[k - 1] - 1))
      
    }
    
    phiM[l] = prod(prod_vec) * 
      exp(lambda * (1 - A[min(which(A[1,] > 0)), 1]) * 
            (h_1 - 1)) *
      lambda * h_1 # seule différence vs algo 3
    
  }
  
  # output pdf_M
  Re(fft(phiM, inverse = T))/nfft
  
}
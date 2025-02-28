###
### Travail 2, ACT-7119
### Theorem 7 : Distribution of the Sum
###

pdf_M = function(A, lambda){
  # A : matrice d'adjacence
  # lambda : paramètre des lois de Poisson
  
  d = nrow(A)
  
  nfft = 2^15
  
  b = c(0, 1, rep(0, nfft - 2)) # fmp v.a. dégénérée à 1
  
  phib = fft(b)
  
  phiM = numeric(nfft)
  
  k_domain = rev(seq(d)[-1])
  
  for (l in seq(nfft)){
    
    H = matrix(data = 1,
               nrow = d,
               ncol = d) # all-1 matrix
    
    prod_vec = numeric(length(k_domain))
    
    for (k in k_domain){
      
      pi_k = min(which(A[k,] > 0)) # infimum (devrait jamais donner infini parce
      # que sinon ça veut dire qu'il y a un node indépendant de tous les autres)
      
      h_k = phib[l] * prod(H[k,])
      
      H[pi_k, k] = 1 - A[pi_k, k] + A[pi_k, k] * h_k # overwrite H
      
      prod_vec[k - 1] = exp(lambda * (1 - A[pi_k, k]) * 
                              (h_k - 1))
      
    }
    
    h_1 = phib[l] * prod(H[1,])
    
    phiM[l] = prod(prod_vec) * exp(lambda * (h_1 - 1)) # ajouter k = 1
    
  }
  
  # output pdf_M
  Re(fft(phiM, inverse = T))/nfft
  
}

# Validation pdf_M ----
alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances

# Matrice adjacente
A = matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

lambda = 1 # paramètre des lois de Poisson

fM = pdf_M(A, lambda)
sum(fM) # somme à 1

m = seq(nfft) - 1
esp_M = sum(m*fM) # changer les dépendances ne devrait pas avoir d'impact
var_M = sum(m^2*fM) - esp_M^2 # augmenter la dépendance devrait augmenter la variance


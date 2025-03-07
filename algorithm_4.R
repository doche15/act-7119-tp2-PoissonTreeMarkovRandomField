###
### Travail 2, ACT-7119
### Algorithm 5 : Computing expected allocations of Nv to M.
###

source("rerooting.R")
source("sample_MPMRF.R")

exp_alloc_Nv = function(A, root_node, lambda, v, nfft = 2^15){
  # A : matrice d'adjacence
  # root_node : numéro de la racine initiale de A
  # lambda : paramètre des lois de Poisson
  # v : numéro du noeud dont on veut l'allocation

  Aprime <-  reroot(A, v)

  d = nrow(Aprime)

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

          pi_k = min(which(Aprime[k,] > 0)) # infimum (devrait jamais donner infini parce
          # que sinon ça veut dire qu'il y a un node indépendant de tous les autres)

          h_k = phib[l] * prod(H[k,])

          H[pi_k, k] = 1 - Aprime[pi_k, k] + Aprime[pi_k, k] * h_k # overwrite H

          prod_vec[k - 1] = exp(lambda * (1 - Aprime[pi_k, k]) *
                                    (h_k - 1))

      }

      h_1 = phib[l] * prod(H[1,])

      phiM[l] = prod(prod_vec) * exp(lambda * (h_1 - 1)) * lambda * h_1 # diff avec algo3

  }

  # output pdf_M
  Re(fft(phiM, inverse = T))/nfft

}


### Validation
alpha12 <-  0.2 ; alpha23 <-  0.4; alpha24 <-  0.7 # dépendances
lambda = 4

# Matrice adjacente
A <-  matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

rea <- rMPMRF(10000000, A, lambda = lambda)

reaM <- rowSums(rea)

expec1k <- exp_alloc_Nv(A, 1, lambda, 4)
m <- 22
expec1k[m]

mean(rea[, 4] * (reaM == (m-1))) # semble OK

###
### Travail 2, ACT-7119
### Classe du modèle MPMRF de [Côté et al., 2025]
###
##


import numpy as np


# Référence pour les classes Python :
# https://docs.python.org/3/tutorial/classes.html


class MPMRF():
    def __init__(self, p_weighted_adjacency_matrix: np.ndarray, p_lambda):
        self.weighted_adjacency_matrix = p_weighted_adjacency_matrix
        self.lambdap = p_lambda
        self.dimension = self.weighted_adjacency_matrix.shape[0]
    
    def sample(self, p_nsamples, p_seed = 11):
        np.random.default_rng(seed=p_seed)
        samples = np.zeros((p_nsamples, self.dimension))

        for i in range(p_nsamples):
            samples[i, 0] = np.random.poisson(self.lambdap)

            for k in range(1, self.dimension):
                pi_k = np.nonzero(self.weighted_adjacency_matrix[k,:])[0][0]

                B_k = np.random.binomial(samples[i, k - 1], self.weighted_adjacency_matrix[pi_k, k])
                
                # Valider ce qu'est le a_k dans l'Algorithme 2 de [Côté et al., 2025]
                L_k = np.random.poisson(self.lambdap * (1 - self.weighted_adjacency_matrix[pi_k, k]))

                samples[i, k] = B_k + L_k
        
        return samples


    def pmf(self, x):
        pass
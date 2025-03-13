#### ACT-7119  - Travail pratique 2
### Tree-structured Markov Random Fields with Poisson Marginal Distribution
###

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

### Reproduction de l'exemple numérique (Section 8)

lambda <- 1

#### Construction de la matrice d'adjacence
beta <- c(0.3, 0.7, 0.9)

d <- 50

## Arètes (selon figure 5 de l'article)
# (1, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8),
# (2, 9), (9, 10), (9, 11), (9, 12), (9, 13), (9, 14), (9, 15),
# (9, 16), (16, 17), (16, 18), (16, 19), (16, 20), (16, 21), (16, 22),
# (16, 23), (23, 24), (23, 25), (23, 26), (23, 27), (23, 28), (23, 29),
# (23, 30), (30, 31), (30, 32), (30, 33), (30, 34), (30, 35), (30, 36),
# (30, 37), (30, 38), (30, 39), (30, 40), (30, 41), (30, 42), (30, 43),
# (30, 44), (30, 45), (30, 46), (30, 47), (30, 48), (30, 49), (30, 50)
d <- 251
A <- diag(d)
lolipop <- c(1, 2, 2, 3, 2, 4, 2, 5, 2, 6, 2, 7, 2, 8,
             2, 9, 9, 10, 9, 11, 9, 12, 9, 13, 9, 14, 9, 15,
             9, 16, 16, 17, 16, 18, 16, 19, 16, 20, 16, 21, 16, 22,
             16, 23, 23, 24, 23, 25, 23, 26, 23, 27, 23, 28, 23, 29,
             23, 30, 30, 31, 30, 32, 30, 33, 30, 34, 30, 35, 30, 36,
             30, 37, 30, 38, 30, 39, 30, 40, 30, 41, 30, 42, 30, 43,
             30, 44, 30, 45, 30, 46, 30, 47, 30, 48, 30, 49, 30, 50)

aretes <- matrix(c(1, 2, lolipop + 1, 1, 52, lolipop + 51,
                   1, 102, lolipop + 101, 1, 152, lolipop + 151,
                   1, 202, lolipop + 201),
                     byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, 0.5)


for (i in 1:nrow(aretes))
    {
        u <- aretes_alpha[i, 1]
        v <- aretes_alpha[i, 2]

        A[u, v] <- aretes_alpha[i, 3]
        A[v, u] <- aretes_alpha[i, 3]
    }



source("graph_A.R")

see_tree_graph(A0.7, 1, titre = "Figure 5")




### pdf_M (algorithme 3)
source("pdf_M.R")
nfft <- 2^15
fM0 <- dpois(0:(nfft - 1), lambda * d)
fM0.3 <- pdf_M(A0.3, lambda)
fM0.7 <- pdf_M(A0.7, lambda)
fM0.9 <- pdf_M(A0.9, lambda)


df_FM <- data.frame(x = 0:(nfft - 1),
                    "beta0" = fM0,
                    "beta0.3" = fM0.3,
                    "beta0.7" = fM0.7,
                    "beta0.9" = fM0.9) %>% pivot_longer(cols = -x,
                                                        names_to = "beta",
                                                        values_to = "fM")


dark2 <- brewer.pal(8, "Dark2")
couleurs <- c("M0" = dark2[1], "M0.3" = dark2[2],
              "M0.7" = dark2[3],
              "M0.9" = dark2[8])

(ggFM <- ggplot(df_FM %>% filter(x <= 200), aes(x = x, y = fM, col = beta)) +
    geom_segment(aes(xend = x, yend = 0), show.legend = TRUE) +
    geom_point() +
    labs(title = bquote("FMP de M selon les valeurs de" ~ beta),
         x = expression(x),
         y = expression(p[M](x)),
         col = expression(beta)) +
    theme_classic() +
        scale_color_manual(values = setNames(couleurs, c("beta0",
                                                         "beta0.3",
                                                         "beta0.7",
                                                         "beta0.9")),
                           labels = c("beta0" = 0,
                                                         "beta0.3" = 0.3,
                                                         "beta0.7" = 0.7,
                                                         "beta0.9" = 0.9)))
ggsave("illustrations_numeriques/figs/repro_section8/fM.pdf",
       width = 9,
       height = 4.8)
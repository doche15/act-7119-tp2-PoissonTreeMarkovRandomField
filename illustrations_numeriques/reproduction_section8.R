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


make_adjmatrix <- function(d, beta)
{
    A <- diag(d)
    aretes <- matrix(c(1, 2, 2, 3, 2, 4, 2, 5, 2, 6, 2, 7, 2, 8,
                       2, 9, 9, 10, 9, 11, 9, 12, 9, 13, 9, 14, 9, 15,
                       9, 16, 16, 17, 16, 18, 16, 19, 16, 20, 16, 21, 16, 22,
                       16, 23, 23, 24, 23, 25, 23, 26, 23, 27, 23, 28, 23, 29,
                       23, 30, 30, 31, 30, 32, 30, 33, 30, 34, 30, 35, 30, 36,
                       30, 37, 30, 38, 30, 39, 30, 40, 30, 41, 30, 42, 30, 43,
                       30, 44, 30, 45, 30, 46, 30, 47, 30, 48, 30, 49, 30, 50),
                     byrow = TRUE, ncol = 2)

    aretes_alpha <- cbind(aretes, beta)


    for (i in 1:nrow(aretes))
    {
        u <- aretes_alpha[i, 1]
        v <- aretes_alpha[i, 2]

        A[u, v] <- aretes_alpha[i, 3]
        A[v, u] <- aretes_alpha[i, 3]
    }
    A
}


A0.3 <- make_adjmatrix(d, beta[1])
A0.7 <- make_adjmatrix(d, beta[2])
A0.9 <- make_adjmatrix(d, beta[3])

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


### Secondary distribution C_M
# ** attention long à rouler, doncje les laisse commentés
source("pgf_Cm.R")
# nfft <- 2^2
# fc0 <- c(0, 1, rep(0, nfft -2))
# ffc0 <- fft(fc0)
# fC0.3 <- Re(fft(sapply(1:nfft, function(a) pgf_Cm(A0.3, lambda, ffc0[a], 1)),
#                 inverse = TRUE))/nfft
#
# fC0.7 <- Re(fft(sapply(1:nfft, function(a) pgf_Cm(A0.7, lambda, ffc0[a], 1)),
#                 inverse = TRUE))/nfft
#
# fC0.9 <- Re(fft(sapply(1:nfft, function(a) pgf_Cm(A0.9, lambda, ffc0[a], 1)),
#                 inverse = TRUE))/nfft
#
# sum(fC0.7)
# sum(fC0.9)
#
# df_FC <- data.frame(x = 0:(nfft - 1),
#                     "beta0" = fc0,
#                     "beta0.3" = fC0.3,
#                     "beta0.7" = fC0.7,
#                     "beta0.9" = fC0.9) %>% pivot_longer(cols = -x,
#                                                         names_to = "beta",
#                                                         values_to = "fM")
#
#
# dark2 <- brewer.pal(8, "Dark2")
# couleurs <- c("M0" = dark2[1], "M0.3" = dark2[2],
#               "M0.7" = dark2[3],
#               "M0.9" = dark2[8])
#
# (ggFM <- ggplot(df_FC %>% filter(x <= 200), aes(x = x, y = fM, col = beta)) +
#         geom_segment(aes(xend = x, yend = 0), show.legend = TRUE) +
#         geom_point() +
#         labs(title = bquote("FMP de M selon les valeurs de" ~ beta),
#              x = expression(x),
#              y = expression(p[M](x)),
#              col = expression(beta)) +
#         theme_classic() +
#         scale_color_manual(values = setNames(couleurs, c("beta0",
#                                                          "beta0.3",
#                                                          "beta0.7",
#                                                          "beta0.9")),
#                            labels = c("beta0" = 0,
#                                       "beta0.3" = 0.3,
#                                       "beta0.7" = 0.7,
#                                       "beta0.9" = 0.9)))



aretes <- matrix(c(1, 2, 2, 3, 2, 4, 2, 5, 2, 6, 2, 7, 2, 8,
                   2, 9, 9, 10, 9, 11, 9, 12, 9, 13, 9, 14, 9, 15,
                   9, 16, 16, 17, 16, 18, 16, 19, 16, 20, 16, 21, 16, 22,
                   16, 23, 23, 24, 23, 25, 23, 26, 23, 27, 23, 28, 23, 29,
                   23, 30, 30, 31, 30, 32, 30, 33, 30, 34, 30, 35, 30, 36,
                   30, 37, 30, 38, 30, 39, 30, 40, 30, 41, 30, 42, 30, 43,
                   30, 44, 30, 45, 30, 46, 30, 47, 30, 48, 30, 49, 30, 50),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, 0)
dd <- nrow(aretes_alpha)

lamM <- function(beta, d, dd, lambda)
{
    lambda * (d - dd * beta)
}

ecm <- function(beta, d, dd)
{
    d / (d - dd * beta)
}

sapply(beta, lamM, d = d, dd = dd, lambda = lambda)
sapply(beta, ecm, d = d, dd = dd)

### Mesures de risque sur fM
nfft <- 2^15
domm <- 0:(nfft - 1)

(EM <- c(sum(domm * fM0),
        sum(domm * fM0.3),
        sum(domm * fM0.7),
        sum(domm * fM0.9)))

(VM <- c(sum(domm^2 * fM0) - sum(domm * fM0)^2,
        sum(domm^2 * fM0.3) - sum(domm * fM0.3)^2,
        sum(domm^2 * fM0.7) - sum(domm * fM0.7)^2,
        sum(domm^2 * fM0.9) - sum(domm * fM0.9)^2))

kappa <- c(0.9, 0.99)

(v0.9 <- c(min(which(cumsum(fM0) >= kappa[1])) - 1,
           min(which(cumsum(fM0.3) >= kappa[1])) - 1,
           min(which(cumsum(fM0.7) >= kappa[1])) - 1,
           min(which(cumsum(fM0.9) >= kappa[1])) - 1))

(tv0.9 <- c(sum(fM0 * pmax(domm - v0.9[1], 0)) / (1 - kappa[1]) + v0.9[1],
            sum(fM0.3 * pmax(domm - v0.9[2], 0)) / (1 - kappa[1]) + v0.9[2],
            sum(fM0.7 * pmax(domm - v0.9[3], 0)) / (1 - kappa[1]) + v0.9[3],
            sum(fM0.9 * pmax(domm - v0.9[4], 0)) / (1 - kappa[1]) + v0.9[4]))

(v0.99 <- c(min(which(cumsum(fM0) >= kappa[2])) - 1,
           min(which(cumsum(fM0.3) >= kappa[2])) - 1,
           min(which(cumsum(fM0.7) >= kappa[2])) - 1,
           min(which(cumsum(fM0.9) >= kappa[2])) - 1))

(tv0.99 <- c(sum(fM0 * pmax(domm - v0.99[1], 0)) / (1 - kappa[2]) + v0.99[1],
            sum(fM0.3 * pmax(domm - v0.99[2], 0)) / (1 - kappa[2]) + v0.99[2],
            sum(fM0.7 * pmax(domm - v0.99[3], 0)) / (1 - kappa[2]) + v0.99[3],
            sum(fM0.9 * pmax(domm - v0.99[4], 0)) / (1 - kappa[2]) + v0.99[4]))

rho <- 0.1

data.frame("beta" = c(0, beta),
           "VarM" = VM,
           "VaR0.9" = v0.9,
           "TVaR0.9" = tv0.9,
           "VaR0.99" = v0.99,
           "TVaR0.99" = tv0.99) # Tableau 2

### Graphiques des fonctions de répartition

df_FM <- data.frame(x = 0:(nfft - 1),
                    "beta0" = cumsum(fM0),
                    "beta0.3" = cumsum(fM0.3),
                    "beta0.7" = cumsum(fM0.7),
                    "beta0.9" = cumsum(fM0.9)) %>% pivot_longer(cols = -x,
                                                        names_to = "beta",
                                                        values_to = "FM")


(ggFM <- ggplot(df_FM %>% filter(x <= 200), aes(x = x, y = FM, col = beta)) +
        geom_line() +
        labs(title = bquote("Fonction de répartition de M selon les valeurs de" ~ beta),
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
ggsave("illustrations_numeriques/figs/repro_section8/frepM.pdf",
       width = 9,
       height = 4.8)



### Contributions à la TVaR

# equation de l'article
contrib_tvar <- function(fM, alloc, lambda, kappa)
{
    FM <- cumsum(fM)
    varr <- min(which(FM >= kappa)) - 1
    aa <- lambda - sum(alloc[0:varr + 1]) + (FM[varr + 1] - kappa) /fM[varr + 1] *
        alloc[varr + 1]

    1 / (1 - kappa) * aa
}
source("algorithm_4.R")
vv <- c(1, 16, 30)

contrib_0 <- tv0.9[1] / d
alloc1_0.3 <- exp_alloc_Nv(A0.3, 1, lambda, 1)
alloc1_0.7 <- exp_alloc_Nv(A0.7, 1, lambda, 1)
alloc1_0.9 <- exp_alloc_Nv(A0.9, 1, lambda, 1)

contrib_0 <- tv0.9[1] / d
alloc16_0.3 <- exp_alloc_Nv(A0.3, 1, lambda, 16)
alloc16_0.7 <- exp_alloc_Nv(A0.7, 1, lambda, 16)
alloc16_0.9 <- exp_alloc_Nv(A0.9, 1, lambda, 16)

contrib_0 <- tv0.9[1] / d
alloc30_0.3 <- exp_alloc_Nv(A0.3, 1, lambda, 30)
alloc30_0.7 <- exp_alloc_Nv(A0.7, 1, lambda, 30)
alloc30_0.9 <- exp_alloc_Nv(A0.9, 1, lambda, 30)

(contrib1 <- c(contrib_0,
              contrib_tvar(fM0.3, alloc1_0.3, lambda, kappa[1]),
              contrib_tvar(fM0.7, alloc1_0.7, lambda, kappa[1]),
              contrib_tvar(fM0.9, alloc1_0.9, lambda, kappa[1])))


(contrib16 <- c(contrib_0,
               contrib_tvar(fM0.3, alloc16_0.3, lambda, kappa[1]),
               contrib_tvar(fM0.7, alloc16_0.7, lambda, kappa[1]),
               contrib_tvar(fM0.9, alloc16_0.9, lambda, kappa[1])))

(contrib30 <- c(contrib_0,
               contrib_tvar(fM0.3, alloc30_0.3, lambda, kappa[1]),
               contrib_tvar(fM0.7, alloc30_0.7, lambda, kappa[1]),
               contrib_tvar(fM0.9, alloc30_0.9, lambda, kappa[1])))

data.frame("beta" = c(0, beta),
           "C(N1, M)" = contrib1,
           "C(N16, M)" = contrib16,
           "C(N30, M)" = contrib30) # Tableau 3

### Changement de la dépendance
lambda <- 0.1
aretes <- matrix(c(1, 2, 2, 3, 2, 4, 2, 5, 2, 6, 2, 7, 2, 8,
                   2, 9, 9, 10, 9, 11, 9, 12, 9, 13, 9, 14, 9, 15,
                   9, 16, 16, 17, 16, 18, 16, 19, 16, 20, 16, 21, 16, 22,
                   16, 23, 23, 24, 23, 25, 23, 26, 23, 27, 23, 28, 23, 29,
                   23, 30, 30, 31, 30, 32, 30, 33, 30, 34, 30, 35, 30, 36,
                   30, 37, 30, 38, 30, 39, 30, 40, 30, 41, 30, 42, 30, 43,
                   30, 44, 30, 45, 30, 46, 30, 47, 30, 48, 30, 49, 30, 50),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, c(rep(0.8, 7),
                                0.1,
                                rep(0.6, 6),
                                0.1,
                                rep(0.5, 6),
                                0.1,
                                rep(0.4, 6),
                                0.1,
                                rep(0.9, 20))
)

A <- diag(d)

for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A[u, v] <- aretes_alpha[i, 3]
    A[v, u] <- aretes_alpha[i, 3]
}



fM <- pdf_M(A, lambda)

df_FM <- data.frame(x = 0:(nfft - 1),
                    "fM" = fM) %>% pivot_longer(cols = -x,
                                                        names_to = "beta",
                                                        values_to = "fM")


dark2 <- brewer.pal(8, "Dark2")
couleurs <- c("fM" = dark2[1])

(ggFM <- ggplot(df_FM %>% filter(x <= 40), aes(x = x, y = fM)) +
        geom_segment(aes(xend = x, yend = 0), show.legend = TRUE,
                     col = couleurs[1]) +
        geom_point(col = couleurs[1]) +
        labs(title = "FMP de M",
             x = expression(x),
             y = expression(p[M](x)),
             col = expression(beta)) +
        theme_classic()) #  Figure 9


### Allocations

alloc1 <- exp_alloc_Nv(A, 1, lambda, 1)
alloc16 <- exp_alloc_Nv(A, 1, lambda, 16)
alloc30 <- exp_alloc_Nv(A, 1, lambda, 30)

df_alloc <- data.frame(x = 0:(nfft - 1),
                    "alloc1" = alloc1 / fM,
                    "alloc16" = alloc16 / fM,
                    "alloc30" = alloc30 / fM) %>% pivot_longer(cols = -x,
                                                names_to = "nu",
                                                values_to = "alloc")


dark2 <- brewer.pal(8, "Dark2")
couleurs <- c("alloc1" = dark2[1],
              "alloc16" = dark2[2],
              "alloc30" = dark2[3])

(ggFM <- ggplot(df_alloc %>% filter(x <= 50), aes(x = x, y = alloc, col = nu)) +
        geom_segment(aes(xend = x, yend = 0), show.legend = TRUE) +
        geom_point() +
        labs(title = bquote("Espérance conditionnelle de" ~ N[nu] ~ "sachant" ~  M == k),
             x = expression(k),
             y = bquote(E[N[nu] ~ "|" ~ M == k]),
             col = expression(nu)) +
        theme_classic() +
    scale_color_manual(values = setNames(couleurs, c("alloc1",
                                                     "alloc16",
                                                     "alloc30")),
                       labels = c("alloc1" = 1,
                                  "alloc16" = 16,
                                  "alloc30" = 30)))#  Figure 10
ggsave("illustrations_numeriques/figs/repro_section8/ENnuMk.pdf",
       width = 9,
       height = 4.8)

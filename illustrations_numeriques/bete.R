#### ACT-7119  - Travail pratique 2
### Tree-structured Markov Random Fields with Poisson Marginal Distribution
###

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

### Reproduction de l'exemple numérique (Section 8)

lambda <- 3

#### Construction de la matrice d'adjacence

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

see_tree_graph(A, 1, titre = "Figure 5")




### Résultats
source("pdf_M.R")
source("Algorithm_4.R")
source("sample_MPMRF.R")
nfft <- 2^15
fM <- pdf_M(A, lambda)
mm <- 0:(nfft - 1)
nn <- 5e5
set.seed(80934)
realisations <- rMPMRF(nn, A, lambda)
reaM <- rowSums(realisations)

sum(mm * fM)
sum(mm^2 * fM) - sum(mm * fM)^2
mean(reaM)
mean(reaM^2) - mean(reaM)^2

cumsum(fM)[750]
mean(reaM <= 749)


(v0.9 <- min(which(cumsum(fM) >= 0.9)) - 1)
sort(reaM)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM) / (1 - 0.9) + v0.9
mean(sort(reaM)[(0.9 * nn):nn])

alloc1 <- exp_alloc_Nv(A, 1, lambda, 1, nfft)
alloc2 <- exp_alloc_Nv(A, 1, lambda, 2, nfft)
alloc51 <- exp_alloc_Nv(A, 1, lambda, 51, nfft)

contrib_tvar <- function(fM, alloc, lambda, kappa)
{
    FM <- cumsum(fM)
    varr <- min(which(FM >= kappa)) - 1
    aa <- lambda - sum(alloc[0:varr + 1]) + (FM[varr + 1] - kappa) /fM[varr + 1] *
        alloc[varr + 1]

    1 / (1 - kappa) * aa
}

contrib_tvar(fM, alloc1, lambda, 0.9)
contrib_tvar(fM, alloc2, lambda, 0.9)
contrib_tvar(fM, alloc51, lambda, 0.9)

mean(realisations[, 1] * (reaM > 855)) / (1 - 0.9)
mean(realisations[, 2] * (reaM > 855)) / (1 - 0.9)
mean(realisations[, 51] * (reaM > 855)) / (1 - 0.9)

df <- data.frame(x = 0:(nfft - 1),
                    "fM" = fM)


dark2 <- brewer.pal(8, "Dark2")

(ggFM <- ggplot(df %>% filter(x <= 1300), aes(x = x, y = fM)) +
    #geom_segment(aes(xend = x, yend = 0), show.legend = TRUE, col = dark2[1]) +
    geom_vline(xintercept = 855, color = dark2[2], linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = 896.88, color = dark2[3], linetype = "dashed", linewidth = 1) +
        annotate("text", x = 730, y = 0.006, label = expression(VaR[0.9](M)), parse = TRUE, hjust = 0) +
        annotate("text", x = 910, y = 0.006, label = expression(TVaR[0.9](M)), parse = TRUE, hjust = 0) +
    geom_point(col = dark2[1]) +
    labs(title = bquote("FMP de M pour la bête"),
         x = expression(x),
         y = expression(p[M](x)),
         col = expression(beta)) +
    theme_classic())

ggsave("illustrations_numeriques/figs/bete/fMbete.pdf",
       width = 9,
       height = 4.8)


df_alloc <- data.frame(x = 0:(nfft - 1),
                       "alloc1" = alloc1 / fM,
                       "alloc2" = alloc2 / fM,
                       "alloc51" = alloc51 / fM) %>% pivot_longer(cols = -x,
                                                                  names_to = "nu",
                                                                  values_to = "alloc")

dark2 <- brewer.pal(8, "Dark2")
couleurs <- c("alloc1" = dark2[1],
              "alloc2" = dark2[2],
              "alloc51" = dark2[3])

(ggFM <- ggplot(df_alloc %>% filter(x <= 1500, x > 350), aes(x = x, y = alloc, col = nu)) +
        #geom_segment(aes(xend = x, yend = 0), show.legend = TRUE) +
        geom_point() +
        labs(title = bquote("Espérance conditionnelle de" ~ N[nu] ~ "sachant" ~  M == k),
             x = expression(k),
             y = bquote(E[N[nu] ~ "|" ~ M == k]),
             col = expression(nu)) +
        theme_classic() +
        scale_color_manual(values = setNames(couleurs, c("alloc1",
                                                         "alloc2",
                                                         "alloc51")),
                           labels = c("alloc1" = 1,
                                      "alloc2" = 2,
                                      "alloc51" = 51)))

ggsave("illustrations_numeriques/figs/bete/enmbete.pdf",
       width = 9,
       height = 4.8)

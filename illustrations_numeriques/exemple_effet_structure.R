###
### Travail 2, ACT-7119
### Exemple effet de la structure
###

library(igraph)

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

d <- 25
A <- diag(25)

lam <- 2
alp <- 0.7

source("sample_MPMRF.R")
source("pdf_M.R")
source("graph_A.R")

#### Structure 1
## Arètes

aretes <- matrix(c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10,
                   10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17,
                   17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24,
                   24, 25),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, alp)
A1 <- A
for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A1[u, v] <- aretes_alpha[i, 3]
    A1[v, u] <- aretes_alpha[i, 3]
}


## Simulations
nn <- 5e6
set.seed(324)
realisations <- rMPMRF(nn, A1, lam)

reaM1 <- rowSums(realisations)
rm(realisations)
## Exact
nfft <- 2^15
fM1 <- pdf_M(A1, lam, nfft)

df <- data.frame(x = 0:(nfft - 1),
           "Str1" = fM1) %>% pivot_longer(cols = -x,
                                       names_to = "structure",
                                       values_to = "fM")

mm <- 0:(nfft - 1)

sum(mm * fM1)
mean(reaM1)

sum(mm^2 * fM1) - sum(mm * fM1)^2
mean(reaM1^2) - mean(reaM1)^2

cumsum(fM1)[20]
mean(reaM1 <= 19)


(v0.9 <- min(which(cumsum(fM1) >= 0.9)) - 1)
sort(reaM1)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM1) / (1 - 0.9) + v0.9
mean(sort(reaM1)[(0.9 * nn):nn])


#### Structure 2
## Arètes

aretes <- matrix(c(1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 8, 1, 9, 1, 10, 1, 11,
                   1, 12, 1, 13, 1, 14, 1, 15, 1, 16, 1, 17, 1, 18, 1, 19, 1, 20,
                   1, 21, 1, 22, 1, 23, 1, 24, 1, 25),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, alp)
A2 <- A
for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A2[u, v] <- aretes_alpha[i, 3]
    A2[v, u] <- aretes_alpha[i, 3]
}


## Simulations
nn <- 5e6
set.seed(3256)
realisations <- rMPMRF(nn, A2, lam)

reaM2 <- rowSums(realisations)
rm(realisations)
## Exact
nfft <- 2^15
fM2 <- pdf_M(A2, lam, nfft)

df <- rbind(df, data.frame(x = 0:(nfft - 1),
                 "Str2" = fM2) %>% pivot_longer(cols = -x,
                                                names_to = "structure",
                                                values_to = "fM"))

mm <- 0:(nfft - 1)

sum(mm * fM2)
mean(reaM2)

sum(mm^2 * fM2) - sum(mm * fM2)^2
mean(reaM2^2) - mean(reaM2)^2

cumsum(fM2)[20]
mean(reaM2 <= 19)


(v0.9 <- min(which(cumsum(fM2) >= 0.9)) - 1)
sort(reaM2)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM2) / (1 - 0.9) + v0.9
mean(sort(reaM2)[(0.9 * nn):nn])


#### Structure 3
## Arètes

aretes <- matrix(c(1, 6, 6, 7, 6, 8, 6, 9, 1, 2, 2, 3, 2, 4, 2, 5, 1, 10, 10,
                   11, 10, 12, 10, 13, 1, 14, 14, 15, 14, 16, 14, 17, 1, 18, 18,
                   19, 18, 20, 18, 21, 1, 22, 22, 23, 22, 24, 22, 25),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, alp)
A3 <- A
for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A3[u, v] <- aretes_alpha[i, 3]
    A3[v, u] <- aretes_alpha[i, 3]
}


## Simulations
nn <- 5e6
set.seed(6366)
realisations <- rMPMRF(nn, A3, lam)

reaM3 <- rowSums(realisations)
rm(realisations)
## Exact
nfft <- 2^15
fM3 <- pdf_M(A3, lam, nfft)
df <- rbind(df, data.frame(x = 0:(nfft - 1),
                           "Str3" = fM3) %>% pivot_longer(cols = -x,
                                                          names_to = "structure",
                                                          values_to = "fM"))

mm <- 0:(nfft - 1)

sum(mm * fM3)
mean(reaM3)

sum(mm^2 * fM3) - sum(mm * fM3)^2
mean(reaM3^2) - mean(reaM3)^2

cumsum(fM3)[20]
mean(reaM3 <= 19)


(v0.9 <- min(which(cumsum(fM3) >= 0.9)) - 1)
sort(reaM3)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM3) / (1 - 0.9) + v0.9
mean(sort(reaM3)[(0.9 * nn):nn])


#### Structure 4
## Arètes

aretes <- matrix(c(1, 2, 2, 7, 2, 6, 2, 5, 2, 4, 2, 3,
                   1, 20, 20, 25, 20, 24, 20, 23, 20, 22, 20, 21,
                   1, 14, 14, 15, 14, 16, 14, 17, 14, 18, 14, 19, 1, 8,
                   8, 9, 8, 10, 8, 11, 8, 12, 8, 13),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, alp)
A4 <- A
for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A4[u, v] <- aretes_alpha[i, 3]
    A4[v, u] <- aretes_alpha[i, 3]
}


## Simulations
nn <- 5e6
set.seed(2526)
realisations <- rMPMRF(nn, A4, lam)

reaM4 <- rowSums(realisations)
rm(realisations)
## Exact
nfft <- 2^15
fM4 <- pdf_M(A4, lam, nfft)
df <- rbind(df, data.frame(x = 0:(nfft - 1),
                           "Str4" = fM4) %>% pivot_longer(cols = -x,
                                                          names_to = "structure",
                                                          values_to = "fM"))

mm <- 0:(nfft - 1)

sum(mm * fM4)
mean(reaM4)

sum(mm^2 * fM4) - sum(mm * fM4)^2
mean(reaM4^2) - mean(reaM4)^2

cumsum(fM4)[20]
mean(reaM4 <= 19)


(v0.9 <- min(which(cumsum(fM4) >= 0.9)) - 1)
sort(reaM4)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM4) / (1 - 0.9) + v0.9
mean(sort(reaM4)[(0.9 * nn):nn])


#### Structure 5
## Arètes

aretes <- matrix(c(1, 2, 2, 3, 2, 4, 4, 5, 5, 6, 5, 7, 4, 8, 4, 9, 4, 10, 4, 11,
                   11, 12, 11, 13, 1, 14, 14, 15, 14, 16, 16, 17, 17, 18, 17, 19,
                   16, 20, 16, 21, 16, 22, 16, 23, 23, 24, 23, 25),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, alp)
A5 <- A
for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A5[u, v] <- aretes_alpha[i, 3]
    A5[v, u] <- aretes_alpha[i, 3]
}


## Simulations
nn <- 5e6
set.seed(2526)
realisations <- rMPMRF(nn, A5, lam)

reaM5 <- rowSums(realisations)
rm(realisations)
## Exact
nfft <- 2^15
df <- rbind(df, data.frame(x = 0:(nfft - 1),
                           "Str5" = fM5) %>% pivot_longer(cols = -x,
                                                          names_to = "structure",
                                                          values_to = "fM"))

mm <- 0:(nfft - 1)

sum(mm * fM5)
mean(reaM5)

sum(mm^2 * fM5) - sum(mm * fM5)^2
mean(reaM5^2) - mean(reaM5)^2

cumsum(fM5)[20]
mean(reaM5 <= 19)


(v0.9 <- min(which(cumsum(fM5) >= 0.9)) - 1)
sort(reaM5)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM5) / (1 - 0.9) + v0.9
mean(sort(reaM5)[(0.9 * nn):nn])


#### Structure 6
## Arètes

aretes <- matrix(c(1, 2, 2, 3, 3, 4, 3, 5, 2, 6, 6, 7, 6, 8, 1, 9, 9, 10,
                   10, 11, 10, 12, 9, 23, 23, 24, 23, 25, 9, 13, 13, 14, 13, 15,
                   1, 16, 16, 17, 17, 18, 17, 19, 16, 20, 20, 21, 20, 22),
                 byrow = TRUE, ncol = 2)

aretes_alpha <- cbind(aretes, alp)
A6 <- A
for (i in 1:nrow(aretes))
{
    u <- aretes_alpha[i, 1]
    v <- aretes_alpha[i, 2]

    A6[u, v] <- aretes_alpha[i, 3]
    A6[v, u] <- aretes_alpha[i, 3]
}


## Simulations
nn <- 5e6
set.seed(22452526)
realisations <- rMPMRF(nn, A6, lam)

reaM6 <- rowSums(realisations)
rm(realisations)
## Exact
nfft <- 2^15
df <- rbind(df, data.frame(x = 0:(nfft - 1),
                           "Str6" = fM6) %>% pivot_longer(cols = -x,
                                                          names_to = "structure",
                                                          values_to = "fM"))

mm <- 0:(nfft - 1)

sum(mm * fM6)
mean(reaM6)

sum(mm^2 * fM6) - sum(mm * fM6)^2
mean(reaM6^2) - mean(reaM6)^2

cumsum(fM6)[20]
mean(reaM6 <= 19)


(v0.9 <- min(which(cumsum(fM6) >= 0.9)) - 1)
sort(reaM6)[0.9 * nn]

sum(pmax(mm - v0.9, 0) * fM6) / (1 - 0.9) + v0.9
mean(sort(reaM6)[(0.9 * nn):nn])


### ggplot

dark2 <- brewer.pal(8, "Dark2")
couleurs <- c("fM" = dark2[1])

(ggFM <- ggplot(df %>% filter(x <= 140, structure %in% c("Str1", "Str2", "Str3", "Str4")), aes(x = x, y = fM, col = structure)) +
    geom_segment(aes(xend = x, yend = 0), show.legend = TRUE) +
    geom_line() +
    labs(title = bquote("Fonction de répartition de M selon la structure"),
         x = expression(x),
         y = expression(p[M](x)),
         col = "Structure") +
    theme_classic() +
    scale_color_manual(values = setNames(dark2[c(8, 3, 6, 4, 2, 1)], c("Str1",
                                                     "Str2",
                                                     "Str3",
                                                     "Str4", "Str5",
                                                     "Str6")),
                       labels = c("Str1" = 1,
                                  "Str2" = 2,
                                  "Str3" = 3,
                                  "Str4" = 4,
                                  "Str5" = 5,
                                  "Str6" = 6)))


ggsave("illustrations_numeriques/figs/effet_structure/fM1234.pdf",
       width = 9,
       height = 4.8)

(ggFM2 <- ggplot(df %>% filter(x <= 140, structure %in% c("Str1", "Str5", "Str6")), aes(x = x, y = fM, col = structure)) +
        geom_segment(aes(xend = x, yend = 0), show.legend = TRUE) +
        geom_line() +
        labs(title = bquote("Fonction de répartition de M selon la structure"),
             x = expression(x),
             y = expression(p[M](x)),
             col = "Structure") +
        theme_classic() +
        scale_color_manual(values = setNames(dark2[c(8, 3, 6, 4, 2, 1)], c("Str1",
                                                                           "Str2",
                                                                           "Str3",
                                                                           "Str4", "Str5",
                                                                           "Str6")),
                           labels = c("Str1" = 1,
                                      "Str2" = 2,
                                      "Str3" = 3,
                                      "Str4" = 4,
                                      "Str5" = 5,
                                      "Str6" = 6)))

ggsave("illustrations_numeriques/figs/effet_structure/fM156.pdf",
       width = 9,
       height = 4.8)

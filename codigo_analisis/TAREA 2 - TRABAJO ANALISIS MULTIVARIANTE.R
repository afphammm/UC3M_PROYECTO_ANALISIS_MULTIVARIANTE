# TAREA 2 - TRABAJO ANÁLISIS MULTIVARIANTE
#---------------------------------------------------------------------------------------------------
library(ggplot2)
library(GGally)
library(corrplot)
library(geigen)
load("EPFhogar_2024.RData")

datos <- Microdatos

datos_sub <- datos[, c(
  "GASTOTGR01", "GASTOTGR04", "GASTOTGR06",  "GASTOTGR09",
  "GASTOTGR11", "GASTOT","OTROIN", 
  "CAJENA","ESTUDREDSP", "TAMAMU", "TAMANO"
)]

colnames(datos_sub) <- c(
  "Gasto_Comida", "Gasto_Vivienda", "Gasto_Sanidad", "Gasto_Ocio",
  "Gasto_Restaurante", "Gasto_total", "Otros_ingresos", 
  "Asalariado", "Estudios", "Tam_Municipio", "Tam_Hogar"
)


datos_limpios <- datos_sub

datos_limpios$Otros_bin <- ifelse(datos_limpios$Otros_ingresos == 1, 1, 0)
datos_limpios$Asalariado_bin <- ifelse(datos_limpios$Asalariado == 1, 1, 0)


datos_limpios$Estudios <- as.factor(datos_limpios$Estudios)
datos_limpios$Tam_Municipio <- as.factor(datos_limpios$Tam_Municipio)
datos_limpios$Tam_Hogar <- as.factor(datos_limpios$Tam_Hogar)

datos_limpios$Asalariado <- NULL
datos_limpios$Otros_ingresos <- NULL

set.seed(1234)
filas <- sample(nrow(datos_limpios), 1000)
datos_final <- datos_limpios[filas, ]

str(datos_final)
summary(datos_final)

datos_num <- datos_final[, sapply(datos_final, is.numeric)]
vars_continuas <- c("Gasto_Comida", "Gasto_Vivienda", 'Gasto_Sanidad', "Gasto_Ocio",
                    'Gasto_Restaurante', "Gasto_total")

datos_num[vars_continuas] <- datos_num[vars_continuas] / 1000

datos_trans <- datos_num
datos_trans[vars_continuas] <- log(datos_trans[vars_continuas]+1)
datos_trans_n <- as.matrix(datos_trans[, vars_continuas])

#----------------------------------------------------------------------------------------------------
# ANÁLISIS CANÓNICO DE POBLACIONES - MANOVA

# Variable multicategórica seleccionada: Nivel de estudios del sustentador principal 
gsize = c(sum(datos_final$Estudios == 1), sum(datos_final$Estudios == 2), 
          sum(datos_final$Estudios == 3), sum(datos_final$Estudios == 4))

grupo = c(rep(1,gsize[1]), rep(2,gsize[2]), rep(3,gsize[3]), rep(4,gsize[4]))

# Gráfico de Dispersión matricial condicional
ggpairs(datos_trans_n, columns = 1:6,
        mapping = aes(color = factor(datos_final$Estudios,
                      levels = 1:4,
                      labels = c("1 - Inferior a Secundaria", "2 - Primera etapa Secundaria",
                                "3 - Segunda etapa Secundaria",  "4 - Educación Superior")),
                      fill = factor(datos_final$Estudios,
                            levels = 1:4,
                            labels = c("1 - Inferior a Secundaria", "2 - Primera etapa Secundaria",
                                     "3 - Segunda etapa Secundaria", "4 - Educación Superior"))),
        upper = list(continuous = wrap("points", alpha = 0.7, size = 0.6)),
        lower = list(continuous = wrap("points", alpha = 0.7, size = 0.6)),
        diag  = list(continuous = wrap("densityDiag", alpha = 0.7)),
        title = "Gráfico de dispersión matricial condicionado al Nivel de Estudios del Sustentador Principal",
        legend = 1                            # <-- esto activa la leyenda
) +
  labs(color = "Nivel de estudios", fill = "Nivel de estudios") +
  theme_bw() +
  theme(
    text  = element_text(size = 16), axis.text  = element_text(size = 12),  
    axis.title = element_text(size = 14), strip.text = element_text(size = 14),
    legend.text = element_text(size = 13),legend.title = element_text(size = 14)
  )
ggsave("dispersion_estudios.png", width = 20, height = 15, dpi = 300)


orden <- order(datos_final$Estudios)
X <- datos_trans_n[orden, ] # datos ordenados según el nivel de estudio

n <- nrow(X)
mmX <- colMeans(X)

# Matriz de centrado global
Hn <- diag(n) - matrix(1, n, n) / n

# Covarianza global
S_global <- t(X) %*% Hn %*% X / n

p <- ncol(X)
g <- length(gsize)
n0 <- cumsum(gsize)
n0_start <- c(1, n0[-g] + 1)  

# Medias por grupo 
mX <- matrix(0, nrow = g, ncol = p)
for (i in 1:g) {
  mX[i, ] <- colMeans(X[n0_start[i]:n0[i], , drop = FALSE])
}

# Matriz de Dispersión entre Grupos
B <- matrix(0, p, p)
for (i in 1:g) {
  diff_i <- matrix(mX[i, ] - mmX, ncol = 1)
  B <- B + gsize[i] * diff_i %*% t(diff_i)
}

# Matriz de Dispersión dentro de los Grupos
W <- matrix(0, p, p)
for (i in 1:g) {
  Xi <- X[n0_start[i]:n0[i], , drop = FALSE]
  ni <- gsize[i]
  Hi <- diag(ni) - matrix(1, ni, ni) / ni
  W  <- W + t(Xi) %*% Hi %*% Xi
}

# Matriz de Dispersión Total
T_mat <- B + W
T_mat2 <- t(X) %*% Hn %*% X   # debe coincidir con B + W

# Test 1 — Igualdad de medias (Lambda de Wilks)
lambda <- det(W) / det(W + B)

wilkstof <- function(L, p, a, b) {
  alpha <- a + b - (p + b + 1) / 2
  beta  <- sqrt((p^2 * b^2 - 4) / (p^2 + b^2 - 5))
  gamma <- (p * b - 2) / 4
  m <- p * b
  n_val <- alpha * beta - 2 * gamma
  n_val <- round(n_val)  
  F_val <- ((1 - L^(1/beta)) / (L^(1/beta))) * n_val / m
  return(list(F = F_val, m = m, n = n_val))
}

wilkstof(lambda, 6, 996, 3)
pf(8.688557,2803,18,lower.tail = FALSE)

# Test 2 — Igualdad de covarianzas (Bartlett)
logH1 <- 0
for (i in 1:g) {
  Xi    <- X[n0_start[i]:n0[i], , drop = FALSE]
  Si    <- t(Xi) %*% (diag(gsize[i]) - matrix(1, gsize[i], gsize[i]) / gsize[i]) %*% Xi / gsize[i]
  logH1 <- logH1 + gsize[i] * log(det(Si))
}

chi     <- abs(n * log(det(S_global)) - logH1)
q       <- (g - 1) * p * (p + 1) / 2
p_valor <- 1 - pchisq(chi, q)


    # Matrices de covarianza por grupo
cov_grupos <- lapply(1:g, function(i) {
  cov(X[n0_start[i]:n0[i], , drop = FALSE])
})

    # Matrices de signos por grupo
signos <- lapply(cov_grupos, sign)

S_pooled   <- W / (n - g)
eig_result <- geigen::geigen(B, S_pooled, symmetric = FALSE)
lambda   <- eig_result$values
idx <- order(lambda, decreasing = TRUE)
lambda  <- lambda[idx]
V   <- eig_result$vectors[, idx]
m <- min(g - 1, p)
V <- V[, 1:m, drop = FALSE]
lambda <- lambda[1:m]
VtSV <- t(V) %*% S_pooled %*% V
V    <- V %*% diag(1 / sqrt(diag(VtSV)))

Y  <- X  %*% V   # coordenadas canónicas de los individuos (n x m)
mY <- mX %*% V   # vectores de medias en nuevas coordenadas  (g x m)


percent <- lambda / sum(lambda) * 100
acum    <- cumsum(percent)


data.frame(
  Eje        = paste0("Can", 1:m),
  Autovalor  = round(lambda, 4),
  Variabilidad = round(percent, 2),
  Acumulado  = round(acum, 2)
)

grupo_vec <- rep(1:g, times = gsize)
df_Y <- data.frame(
  Can1  = Y[, 1],
  Can2  = Y[, 2],
  Grupo = factor(grupo_vec,
                 labels = c("1 - Inferior a Secundaria",
                            "2 - Primera etapa Secundaria",
                            "3 - Segunda etapa Secundaria",
                            "4 - Educación Superior"))
)
df_mY <- data.frame(
  Can1  = mY[, 1],
  Can2  = mY[, 2],
  Grupo = factor(1:g),
  Label = as.character(1:g)
)

ggplot(df_Y, aes(x = Can1, y = Can2, color = Grupo)) +
  geom_point(alpha = 0.5) +
  geom_point(data = df_mY, shape = 17, size = 4, color = "black") +
  geom_text(data = df_mY, aes(label = Grupo), vjust = -1, 
            color = "black", fontface = "bold") +
  labs(
    title = paste0("Coordenadas canónicas (", round(acum[2], 1), "% variabilidad acumulada)"),
    x = paste0("1er eje canónico (", round(percent[1], 1), "%)"),
    y = paste0("2º eje canónico  (", round(percent[2], 1), "%)")
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 16), axis.text  = element_text(size = 12),
    axis.title  = element_text(size = 14), strip.text = element_text(size = 14),  
    legend.text = element_text(size = 13), legend.title = element_text(size = 14)
  )

ggsave("coordenadas_indiv_medios.png", width = 15, height = 10, dpi = 150)


# CÁLCULO DE REGIONES DE CONFIANZA
alpha <- 0.05
radios <- sapply(1:g, function(i) {
  ni <- gsize[i]
  sqrt((p * (n - g)) / ((n - g - p + 1) * ni) * qf(1 - alpha, p, n - g - p + 1))
})

data.frame(Grupo = 1:g, n_i = gsize, Radio = round(radios, 4))

# Puntos del círculo
theta <- seq(0, 2 * pi, length.out = 314)

df_circulos <- data.frame(
  x     = c(outer(cos(theta), radios)) + rep(mY[, 1], each = length(theta)),
  y     = c(outer(sin(theta), radios)) + rep(mY[, 2], each = length(theta)),
  Grupo = factor(rep(1:g, each = length(theta)),
                 labels = c("1 - Inferior a Secundaria","2 - Primera etapa Secundaria",
                            "3 - Segunda etapa Secundaria", "4 - Educación Superior"))
)

ggplot() +
  geom_path(data = df_circulos, aes(x = x, y = y, color = Grupo), 
            linewidth = 1) +
  geom_point(data = df_mY, aes(x = Can1, y = Can2, color = factor(Grupo,
                              labels = c("1 - Inferior a Secundaria", "2 - Primera etapa Secundaria",
                                     "3 - Segunda etapa Secundaria", "4 - Educación Superior"))),
             shape = 17, size = 4) +
  geom_text(data = df_mY, aes(x = Can1, y = Can2, label = Label),
            vjust = -1, fontface = "bold", size = 5) +
  labs(
    title = paste0("Regiones de confianza al 95% (", round(acum[2], 1), "% variabilidad acumulada)"),
    x     = paste0("1er eje canónico (", round(percent[1], 1), "%)"),
    y     = paste0("2º eje canónico  (", round(percent[2], 1), "%)"),
    color = "Nivel de estudios"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 16), axis.text = element_text(size = 12),
    axis.title   = element_text(size = 14), legend.text  = element_text(size = 13),
    legend.title = element_text(size = 14)
  )


ggsave("regiones_confianza.png", width = 15, height = 10, dpi = 150)


data.frame(Grupo = 1:g, Can1 = mY[,1], Can2 = mY[,2])
 

#--------------------------------------------------------------------------------------------------

#PCA

Z_pca <- scale(datos_trans_n, center = TRUE, scale = TRUE)
pca <- prcomp(Z_pca, center = FALSE, scale. = FALSE)

autovalores <- pca$sdev^2
var_exp <- 100 * autovalores / sum(autovalores)
var_acum <- cumsum(var_exp)

tabla_pca <- data.frame(
  Componente = paste0("CP", 1:length(autovalores)),
  Autovalor = round(autovalores, 4),
  Varianza = round(var_exp, 2),
  Acumulada = round(var_acum, 2)
)

tabla_cargas <- data.frame(
  Variable = rownames(pca$rotation),
  CP1 = round(-pca$rotation[, 1], 4),
  CP2 = round(-pca$rotation[, 2], 4)
)

print(tabla_pca)
print(tabla_cargas)

# =========================
# GRÁFICA 1: SCREE PLOT
# =========================
plot(var_exp, type = "b", pch = 19,
     xlab = "Componente principal",
     ylab = "Porcentaje de varianza explicada",
     main = "Scree plot del PCA")
abline(h = 100 / ncol(Z_pca), lty = 2)

# =========================
# GRÁFICA 2: BIPLOT
# =========================
# =========================
# BIPLOT LIMPIO
# =========================

# Coordenadas de individuos
scores <- as.data.frame(pca$x[, 1:2])
colnames(scores) <- c("CP1", "CP2")
scores$Estudios <- datos_final$Estudios

# Coordenadas de variables
loadings <- as.data.frame(-pca$rotation[, 1:2])
colnames(loadings) <- c("CP1", "CP2")
loadings$Variable <- rownames(loadings)

# Escalado de flechas para que se vean bien junto a los puntos
mult <- min(
  (max(scores$CP1) - min(scores$CP1)) / (max(loadings$CP1) - min(loadings$CP1)),
  (max(scores$CP2) - min(scores$CP2)) / (max(loadings$CP2) - min(loadings$CP2))
) * 0.35

loadings$CP1 <- loadings$CP1 * mult
loadings$CP2 <- loadings$CP2 * mult

# Gráfico
ggplot() +
  geom_point(data = scores,
             aes(x = CP1, y = CP2, color = Estudios),
             alpha = 0.35, size = 1.2) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = CP1, yend = CP2),
               arrow = arrow(length = unit(0.2, "cm")),
               linewidth = 0.8) +
  geom_text(data = loadings,
            aes(x = CP1, y = CP2, label = Variable),
            vjust = -0.7, size = 4) +
  labs(
    title = "Biplot del PCA",
    x = paste0("CP1 (", round(var_exp[1], 2), "%)"),
    y = paste0("CP2 (", round(var_exp[2], 2), "%)")
  ) +
  coord_equal() +
  theme_bw()

# =========================
# GRÁFICA 3: INDIVIDUOS EN EL PLANO PRINCIPAL
# =========================
scores <- as.data.frame(pca$x[, 1:2])
colnames(scores) <- c("CP1", "CP2")
scores$Estudios <- datos_final$Estudios

graf_individuos <- ggplot(scores, aes(x = CP1, y = CP2, color = Estudios)) +
  geom_point(alpha = 0.6) +
  labs(
    title = "Hogares en el plano principal",
    x = paste0("CP1 (", round(var_exp[1], 2), "%)"),
    y = paste0("CP2 (", round(var_exp[2], 2), "%)")
  ) +
  coord_equal() +
  theme_bw()

print(graf_individuos)

# =========================
# GRÁFICA 4: CARGAS DE LAS VARIABLES
# =========================
loadings <- as.data.frame(-pca$rotation[, 1:2])
colnames(loadings) <- c("CP1", "CP2")
loadings$Variable <- rownames(loadings)

graf_cargas <- ggplot(loadings) +
  geom_segment(aes(x = 0, y = 0, xend = CP1, yend = CP2),
               arrow = arrow(length = grid::unit(0.2, "cm"))) +
  geom_text(aes(x = CP1, y = CP2, label = Variable), vjust = -0.5) +
  labs(
    title = "Cargas de las variables en el plano principal",
    x = "CP1",
    y = "CP2"
  ) +
  coord_equal() +
  theme_bw()

print(graf_cargas)




# TAREA 3 - TRABAJO ANÁLISIS MULTIVARIANTE
#---------------------------------------------------------------------------------------------------
library(ggplot2)
#library(GGally)
#library(corrplot)
#library(geigen)
library(dbrobust)
library(reshape2)

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
datos_trans$Estudios <- datos_final$Estudios
datos_trans$Tam_Municipio <- datos_final$Tam_Municipio
datos_trans$Tam_Hogar <- datos_final$Tam_Hogar

D2 <- dbrobust::robust_distances(datos_trans, p = c(6,2,3), method = "ggower", return_dist = FALSE)

#----------------------------------------------------------------------------------------------------------

# ANÁLISIS DE COORDENADAS PRINCIPALES - MDS
D <- sqrt(D2)
n <- nrow(D2)
H <- diag(n) - (1/n) * matrix(1, nrow = n, ncol = n)
G <- (-1/2) * H %*% D2 %*% H

L <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
m <- min(L)
epsilon <- 1e-6
if (abs(m) > epsilon) {
  D <- D2 - (2 * m * matrix(1, nrow = n, ncol = n)) + (2 * m * diag(n))
  G <- (-1/2) * H %*% D %*% H
}

MDS <- cmdscale(D, k = 2, eig = TRUE)

Y <- MDS$points
colnames(Y) <- c("Dim1", "Dim2")

vaps <- MDS$eig
percent <- (vaps/sum(abs(vaps)))*100
acum <- cumsum(percent)

print(paste0("Varianza total explicada: ", round(percent[1] + percent[2], 2), "%"))

U = eigen(G)$vectors

n <- nrow(as.matrix(D2))
MDS_completo <- cmdscale(sqrt(D2), k = n - 1, eig = TRUE)
Y_total <- MDS_completo$points
dist_mapa <- as.matrix(dist(Y_total))
dist_original <- sqrt(D2)
diferencia <- max(abs(dist_original - dist_mapa))
print(paste("Diferencia con todas las dimensiones:", diferencia))

  # Representación Gráfica


plot(Y[, 1], Y[, 2], 
     col = "blue", pch = 16, cex = 1.5,
     xlab = "Primera coordenada principal", 
     ylab = "Segunda coordenada principal",
     main = paste0("Porcentaje de variab. explicada ", round(acum[2], 2), "%"),
     panel.first = grid())

if (n <= 100) {
  text(Y[, 1], Y[, 2], labels = 1:n, pos = 3, cex = 0.8)
}

df_var <- data.frame(
  Coord = seq_along(percent),
  Indiv = percent,
  Acum  = acum
)


ggplot(df_var[1:20, ], aes(x = Coord)) +
  geom_col(aes(y = Indiv), fill = "steelblue", alpha = 0.7, width = 0.7) +
  geom_line(aes(y = Acum), color = "firebrick", linewidth = 1.2) +
  geom_point(aes(y = Acum), color = "firebrick", size = 3) +
  geom_hline(yintercept = 80, linetype = 2, color = "grey50", linewidth = 0.8) +
  scale_x_continuous(breaks = 1:20) +
  labs(
    title = "Porcentaje de variabilidad explicada y acumulada",
    subtitle = "Primeras 20 coordenadas principales",
    x = "Coordenada principal",
    y = "% de variabilidad"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave("porcentajes_explicados.png",width = 6, height = 5, dpi = 150)

# CORRELACIONES CRUZADAS

p <- ncol(datos_trans)
k <- 2  ; pcuant   <- 6 ; pnominal <- 2 ; pordinal <- p - pcuant - pnominal

  # --- V de Cramer ---
CramerV <- function(x, y, bins = 4) {
  y_cat <- cut(y, breaks = bins, labels = FALSE, include.lowest = TRUE)
  tbl   <- table(x, y_cat)
  chi2  <- suppressWarnings(chisq.test(tbl)$statistic)
  n     <- sum(tbl)
  k_loc <- min(nrow(tbl), ncol(tbl))
  v     <- sqrt(chi2 / (n * (k_loc - 1)))
  return(as.numeric(v))
}

corr_table <- matrix(0, nrow = p, ncol = k)

    # 1. Pearson para cuantitativas
if (pcuant > 0) {
  corr_table[1:pcuant, ] <- cor(
    datos_trans[, 1:pcuant, drop = FALSE],
    Y[, 1:k, drop = FALSE],          
    method = "pearson"
  )
}

    # 2. V de Cramer para nominales
if (pnominal > 0) {
  idx <- (pcuant + 1):(pcuant + pnominal)
  for (col_i in idx) {
    corr_table[col_i, ] <- sapply(  
      1:k,
      function(j) CramerV(datos_trans[, col_i], Y[, j])
    )
  }
}

    # 3. Spearman para ordinales
if (pordinal > 0) {
  ord_cols <- (pcuant + pnominal + 1):p
  datos_ordinales_num <- data.frame(lapply(datos_trans[, ord_cols], function(x) as.numeric(as.factor(x))))
  corr_table[ord_cols, ] <- cor(datos_ordinales_num, Y[, 1:k], method = "spearman", use = "complete.obs")
}

colnames(corr_table) <- paste0("PC", 1:k)   # "PC1", "PC2" (o hasta k)
rownames(corr_table) <- colnames(datos_trans)

    # --- Heatmap ---
my_colors <- c(
  "#0000CC", "#0066FF", "#66CCFF",
  "#FFFFFF",
  "#FF99FF", "#FF33CC", "#800080"
)


df <- melt(corr_table)
colnames(df) <- c("Variable", "PC", "value")

p_heatmap <- ggplot(df, aes(x = PC, y = Variable, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = my_colors,
    limits = c(-1, 1),
    name   = ""
  ) +
  geom_text(aes(label = round(value, 2)), size = 3) +
  labs(title = "Principal Coordinates Heatmap") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_blank()
  )

print(p_heatmap)

ggsave("pheatmap.png",width = 6, height = 5, dpi = 150)




# --- PERFILES: VARIABLES CUANTITATIVAS ---
umbral    <- 0.3
var_names <- colnames(datos_trans)[1:pcuant]

max_corr <- apply(abs(corr_table[1:pcuant, ]), 1, max)
vars_sel <- which(max_corr >= umbral)
cat("Variables seleccionadas:", names(vars_sel))

if (!dir.exists("plots_cuantis")) dir.create("plots_cuantis")

for (j in vars_sel) {
  df <- data.frame(PCo1 = Y[, 1], PCo2 = Y[, 2], valor = datos_trans[, j])
  
  fig <- ggplot(df, aes(x = PCo1, y = PCo2, color = valor)) +
    geom_point(size = 2.5) +
    scale_color_distiller(palette = "RdYlBu", name = var_names[j]) +
    labs(title = paste("Variable cuantitativa:", var_names[j]),
         x = "PCo1", y = "PCo2") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(fig)
  ggsave(file.path("plots_cuantis", paste0(var_names[j], ".png")),
         plot = fig, width = 6, height = 5, dpi = 150)
}

# --- PERFILES: VARIABLES CUALITATIVAS (nominales + ordinales) ---
idx_cuali    <- (pcuant + 1):p              
var_names_cuali <- colnames(datos_trans)[idx_cuali]

max_corr_cuali <- apply(abs(corr_table[idx_cuali, ]), 1, max)
vars_sel_cuali <- which(max_corr_cuali >= umbral)
cat("Variables seleccionadas:", names(vars_sel_cuali))

if (!dir.exists("plots_cualis")) dir.create("plots_cualis")

for (j in vars_sel_cuali) {
  col_i  <- idx_cuali[j]                          
  grupo  <- factor(datos_trans[, col_i])
  df     <- data.frame(PCo1 = Y[, 1], PCo2 = Y[, 2], grupo = grupo)
  
  fig <- ggplot(df, aes(x = PCo1, y = PCo2, color = grupo)) +
    geom_point(size = 2.5) +
    labs(title = paste("Variable cualitativa:", var_names_cuali[j]),
         x = "PCo1", y = "PCo2",
         color = var_names_cuali[j]) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(fig)
  ggsave(file.path("plots_cualis", paste0(var_names_cuali[j], ".png")),
         plot = fig, width = 6, height = 5, dpi = 150)
  
}

#-----------------------------------------------------------------------------------------
# ANÁLISIS DE CLASIFICACIÓN JERÁRQUICA

# Carpeta donde se guardarán las figuras
carpeta_cluster <- "graficas_clasificacion_estilo_companeros"

if (!dir.exists(carpeta_cluster)) {
  dir.create(carpeta_cluster)
}

# ------------------------------------------------------------------------------------------
# 1. Preparación de tablas para las gráficas
# ------------------------------------------------------------------------------------------

# Como D2 contiene distancias al cuadrado, usamos la raíz cuadrada para hclust
D <- as.dist(sqrt(pmax(D2, 0)))

# Clasificación jerárquica aglomerativa mediante Ward
clasif_jer <- hclust(D, method = "ward.D2")

# Corte en 3 grupos
k <- 3
grupos <- cutree(clasif_jer, k = k)

# Base con los grupos añadidos
datos_perfil <- datos_final
datos_perfil[vars_continuas] <- datos_perfil[vars_continuas] / 1000
datos_perfil$Grupo <- as.factor(grupos)

# Tabla de tamaños de grupo
df_tam <- as.data.frame(table(datos_perfil$Grupo))
colnames(df_tam) <- c("Grupo", "Frecuencia")

# Tabla de medias de gasto por grupo
medias_gasto_grupo <- aggregate(
  datos_perfil[, vars_continuas],
  by = list(Grupo = datos_perfil$Grupo),
  FUN = mean
)

# Tabla de proporciones laborales
prop_otros <- prop.table(table(datos_perfil$Grupo, datos_perfil$Otros_bin), 1)
prop_asal <- prop.table(table(datos_perfil$Grupo, datos_perfil$Asalariado_bin), 1)

df_laboral <- data.frame(
  Grupo = levels(datos_perfil$Grupo),
  Otros_ingresos = prop_otros[, "1"],
  Asalariado = prop_asal[, "1"]
)

# ------------------------------------------------------------------------------------------
# 2. Paleta de colores parecida a la usada en el MDS
# ------------------------------------------------------------------------------------------

col_grupos <- c("#F8766D", "#7CAE00", "#00BFC4")

my_colors <- c(
  "#0000CC", "#0066FF", "#66CCFF",
  "#FFFFFF",
  "#FF99FF", "#FF33CC", "#800080"
)

# ==========================================================================================
# FIGURA 7. Tamaño de los grupos
# ==========================================================================================

p_tam <- ggplot(df_tam, aes(x = Grupo, y = Frecuencia, fill = Grupo)) +
  geom_col(color = "white", width = 0.75) +
  geom_text(aes(label = Frecuencia), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = col_grupos) +
  labs(
    title = "Clasificación jerárquica: tamaño de los grupos",
    x = "Grupo",
    y = "Número de hogares"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

print(p_tam)

ggsave(
  filename = file.path(carpeta_cluster, "07_tamano_grupos_estilo_companeros.png"),
  plot = p_tam,
  width = 6,
  height = 5,
  dpi = 300
)

# ==========================================================================================
# FIGURA 8. Heatmap del perfil medio de gasto por grupo
# ==========================================================================================

# Estandarizamos las medias por variable para comparar perfiles en una misma escala
medias_z <- medias_gasto_grupo

for (v in vars_continuas) {
  medias_z[[v]] <- as.numeric(scale(medias_gasto_grupo[[v]]))
}

# Pasamos a formato largo para ggplot
df_heat_gasto <- melt(
  medias_z,
  id.vars = "Grupo",
  variable.name = "Variable",
  value.name = "Media_estandarizada"
)

p_heat_gasto <- ggplot(
  df_heat_gasto,
  aes(x = Grupo, y = Variable, fill = Media_estandarizada)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Media_estandarizada, 2)), size = 3) +
  scale_fill_gradientn(
    colors = my_colors,
    limits = c(-1.5, 1.5),
    name = "Media\nestandarizada"
  ) +
  labs(title = "Cluster Profile Heatmap") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

print(p_heat_gasto)

ggsave(
  filename = file.path(carpeta_cluster, "08_heatmap_perfil_gasto_grupos.png"),
  plot = p_heat_gasto,
  width = 6,
  height = 6,
  dpi = 300
)

# ==========================================================================================
# FIGURA 9. Perfil medio estandarizado de gasto por grupo
# Esta figura es opcional. Sirve como alternativa al heatmap.
# ==========================================================================================

p_perfil <- ggplot(
  df_heat_gasto,
  aes(
    x = Variable,
    y = Media_estandarizada,
    group = Grupo,
    color = Grupo
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = col_grupos) +
  labs(
    title = "Perfil medio de gasto por grupo",
    x = "Variable de gasto",
    y = "Media estandarizada"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p_perfil)

ggsave(
  filename = file.path(carpeta_cluster, "09_perfil_medio_gasto_estilo_companeros.png"),
  plot = p_perfil,
  width = 8,
  height = 5.5,
  dpi = 300
)

# ==========================================================================================
# FIGURA 10. Heatmap de variables laborales por grupo
# ==========================================================================================

df_laboral_largo <- melt(
  df_laboral,
  id.vars = "Grupo",
  variable.name = "Variable",
  value.name = "Proporcion"
)

# Cambiamos nombres para que salgan más bonitos en la gráfica
df_laboral_largo$Variable <- factor(
  df_laboral_largo$Variable,
  levels = c("Otros_ingresos", "Asalariado"),
  labels = c("Otros ingresos", "Asalariado")
)

p_heat_laboral <- ggplot(
  df_laboral_largo,
  aes(x = Grupo, y = Variable, fill = Proporcion)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Proporcion, 2)), size = 3.5) +
  scale_fill_gradientn(
    colors = my_colors,
    limits = c(0, 1),
    name = "Proporción"
  ) +
  labs(title = "Variables cualitativas por grupo") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

print(p_heat_laboral)

ggsave(
  filename = file.path(carpeta_cluster, "10_heatmap_variables_laborales_grupos.png"),
  plot = p_heat_laboral,
  width = 6,
  height = 4,
  dpi = 300
)


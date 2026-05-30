#----------------------------------------------------------------------------------------------
# CARGAR DE DATOS Y OBTENCIÓN DEL CONJUNTO DE DATOS FINAL
library(GGally)
library(corrplot)
load("EPFhogar_2024.RData")

datos <- Microdatos

  # 2. Selección de las variables clave
datos_sub <- datos[, c(
  "GASTOTGR01", "GASTOTGR04", "GASTOTGR06",  "GASTOTGR09",
  "GASTOTGR11", "GASTOT","OTROIN", 
  "CAJENA","ESTUDIOSSP", "TAMAMU", "TAMANO"
)]

colnames(datos_sub) <- c(
  "Gasto_Comida", "Gasto_Vivienda", "Gasto_Sanidad", "Gasto_Ocio",
  "Gasto_Restaurante", "Gasto_total", "Otros_ingresos", 
  "Asalariado", "Estudios", "Tam_Municipio", "Tam_Hogar"
)

  # 3. Limpieza y Recodificación

datos_limpios <- datos_sub

datos_limpios$Otros_bin <- ifelse(datos_limpios$Otros_ingresos == 1, 1, 0)
datos_limpios$Asalariado_bin <- ifelse(datos_limpios$Asalariado == 1, 1, 0)


datos_limpios$Estudios <- as.factor(datos_limpios$Estudios)
datos_limpios$Tam_Municipio <- as.factor(datos_limpios$Tam_Municipio)
datos_limpios$Tam_Hogar <- as.factor(datos_limpios$Tam_Hogar)

datos_limpios$Asalariado <- NULL
datos_limpios$Otros_ingresos <- NULL

  # 4. Crear el Subset Aleatorio de 1000 filas
set.seed(1234)
filas <- sample(nrow(datos_limpios), 1000)
datos_final <- datos_limpios[filas, ]

str(datos_final)
summary(datos_final)


#----------------------------------------------------------------------------------------------
# EJERCICIO 1 - TAREA 1

datos_num <- datos_final[, sapply(datos_final, is.numeric)]
vars_continuas <- c("Gasto_Comida", "Gasto_Vivienda", 'Gasto_Sanidad', "Gasto_Ocio",
                    'Gasto_Restaurante', "Gasto_total")

datos_num[vars_continuas] <- datos_num[vars_continuas] / 1000
m <- colMeans(datos_num[vars_continuas]) # vector de medias
R <- cor(datos_num[vars_continuas]) # matriz de correlaciones
n <- nrow(datos_num)
S_muestral <- cov(datos_num[vars_continuas])
S <- S_muestral * (n - 1) / n # matriz de covarianzas poblacional

corrplot(R, 
         method = "color",       
         addCoef.col = "black",  
         tl.col = "black",      
         tl.srt = 45,            
         number.cex = 0.8)

ggpairs(datos_num[,vars_continuas], 
        lower = list(continuous = wrap("points", alpha = 1, size = 0.5, color = "blue")), 
        diag = list(continuous = wrap("barDiag", fill = "lightblue")), 
        upper = list(continuous = wrap("points", alpha = 1, size = 0.5, color = "blue"))) +
  theme_minimal() + 
  theme(
    strip.text = element_text(size = 7), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5)
  )  # grafico de dispersión matricial

#ggsave("matriz_dispersion.png", width = 10, height = 10, dpi = 300)


    # Transformaciones no lineales
datos_trans <- datos_num
datos_trans[vars_continuas] <- log(datos_trans[vars_continuas]+1)
datos_trans_n <- as.matrix(datos_trans)

m_new <- colMeans(datos_trans_n[,vars_continuas])
R_new <- cor(datos_trans_n[,vars_continuas])
S_new <- cov(datos_trans_n[,vars_continuas])*(n-1)/n

corrplot(R_new, 
         method = "color",       
         addCoef.col = "black",  
         tl.col = "black",      
         tl.srt = 45,            
         number.cex = 0.8)

ggpairs(datos_trans[vars_continuas], 
        lower = list(continuous = wrap("points", alpha = 1, size = 0.5, color = "blue")), 
        diag = list(continuous = wrap("barDiag", fill = "lightblue")), 
        upper = list(continuous = wrap("points", alpha = 1, size = 0.5, color = "blue"))) +
  theme_minimal() + 
  theme(
    strip.text = element_text(size = 7), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), 
    axis.text.y = element_text(size = 5)
  )  
#ggsave("matriz_normal_dispersion.png", width = 10, height = 10, dpi = 300)

  # Medidas escalares de dispersión y de interdependencia
var_gen <- det(S)
var_gen_new <- det(S_new)
var_total <- sum(diag(S))
var_total_new <- sum(diag(S_new))
coef_interdependencia <- 1-det(R)
coef_interdependencia_new <- 1-det(R_new)


# EJERCICIO 2 - TAREA 1

datos_trans$Asalariado_bin <- as.factor(datos_trans$Asalariado_bin)

ggpairs(
  data = datos_trans,
  columns = vars_continuas,
  mapping = aes(color = Asalariado_bin, fill = Asalariado_bin, alpha = 0.5),
  title = "Matriz de Dispersion Condicional del Gasto por Condicion Laboral",
  upper = list(continuous = wrap("points", alpha = 0.3, size = 0.6)),
  lower = list(continuous = wrap("points", alpha = 0.3, size = 0.6)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.3))
) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 8, face = "bold")
  )

#ggsave("grafico_condicional.png", width = 10, height = 10, dpi = 300)


# 6. Test de Hotelling y supuestos previos
library(ICSNP)

grupo_0 <- datos_trans[datos_trans$Asalariado_bin=="0", vars_continuas]
grupo_1 <- datos_trans[datos_trans$Asalariado_bin=="1", vars_continuas]

test_hotelling <- HotellingsT2(grupo_0, grupo_1)
print(test_hotelling)


# EJERCICIO 3
datos_trans$Estudios <- factor(
  datos_final$Estudios,
  levels = 1:8,
  labels = c(
    "Sin estudios (<5 años)",
    "Educación primaria",
    "Secundaria básica (ESO)",
    "Secundaria superior (Bach/FP medio)",
    "FP grado superior",
    "Universitarios grado corto",
    "Universitarios superiores",
    "Doctorado"
  )
)

ggpairs(
  datos_trans,
  columns = vars_continuas,
  mapping = aes(color = Estudios, fill = Estudios),
  title  = "Análisis de Gasto por Nivel de Estudios",
  legend = 1,  
  upper  = list(continuous = wrap("points", alpha = 0.7, size = 0.6)),
  lower  = list(continuous = wrap("points", alpha = 0.7, size = 0.6)),
  diag   = list(continuous = wrap("densityDiag", alpha = 0.7))
) +
  scale_color_discrete(name = "Nivel de estudios\ndel sustentador principal") +
  scale_fill_discrete(name  = "Nivel de estudios\ndel sustentador principal") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 7),
    legend.title    = element_text(size = 8, face = "bold")
  )
#ggsave("grafico_condicional_estudios.png", width = 10, height = 10, dpi = 300)


# 3b. Contraste de Medias Multivariante: Lambda de Wilks
Y <- as.matrix(datos_trans[, vars_continuas])
G <- datos_trans$Estudios                 

# Ajustamos un modelo lineal multivariante 
fit <- lm(Y ~ G)

# matriz E (Suma de cuadrados del Error / Residuals)
E <- t(residuals(fit)) %*% residuals(fit)

# matriz H (Suma de cuadrados de la Hipótesis / Efecto del grupo)
# matriz Total T = t(Y - MediaTotal) %*% (Y - MediaTotal)
Y_centrado <- scale(Y, scale = FALSE)
T_total <- t(Y_centrado) %*% Y_centrado
H <- T_total - E

# Estadístico Lambda de Wilks (Wilks' Lambda)
# Lambda = |E| / |E + H|
lambda_wilks <- det(E) / det(E + H)

cat("\n--- CONTRASTE DE COMPARACIÓN DE VECTORES DE MEDIAS ---")
cat("\nEstadístico Lambda de Wilks calculado (Manual):", lambda_wilks)

n = 1000
g = 8
p = 6
a = n-g
b = g-1
alpha = a+b-(p+b+1)/2 
beta_2 = (p^2*b^2-4)/(p^2+b^2-5)
gamma = (p*b-2)/4
F_rao = ( (1-lambda_wilks^(1/sqrt(beta_2)))/(lambda_wilks^(1/sqrt(beta_2))) *
            ( (alpha*sqrt(beta_2)-2*gamma)) / (p*b) )

df1 <- p * b                               
df2 <- alpha * sqrt(beta_2) - 2 * gamma       

p_valor <- pf(F_rao, df1 = df1, df2 = df2, lower.tail = FALSE)




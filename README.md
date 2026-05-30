# Análisis Multivariante del Gasto de los Hogares — EPF 2024 (INE)

> Proyecto académico · Grado en Estadística y Empresa · UC3M  
> Herramientas: **R**, **RStudio** | Fuente de datos: **Encuesta de Presupuestos Familiares 2024 — INE**

---

## Descripción

Este proyecto aplica técnicas estadísticas multivariantes sobre los microdatos públicos de la **Encuesta de Presupuestos Familiares (EPF) 2024** del Instituto Nacional de Estadística (INE) para identificar patrones de consumo en los hogares españoles y segmentar perfiles de gasto con potencial aplicación en decisiones de negocio.

---

## Objetivos

- Explorar la estructura del gasto en hogares españoles a partir de datos reales del INE
- Detectar grupos de gasto significativamente distintos mediante análisis estadístico riguroso
- Reducir la dimensionalidad del problema para visualizar patrones latentes
- Segmentar los hogares en perfiles de consumo accionables

---

## Estructura del proyecto

```
epf-analisis-multivariante/
│
├── data/
│   └── raw/             # Microdatos originales del INE (no incluidos por tamaño)
│
├── R/
│   ├── 01_limpieza.R    # Gestión y preprocesado de microdatos
│   ├── 02_eda.R         # Análisis exploratorio y descriptiva
│   ├── 03_hipotesis.R   # Tests de hipótesis multivariantes (MANOVA, Hotelling)
│   ├── 04_pca_mds.R     # Reducción de dimensionalidad (PCA y MDS)
│   └── 05_clustering.R  # Segmentación de hogares (k-means y jerárquico)
│
├── outputs/
│   ├── plots/           # Visualizaciones generadas
│   └── tablas/          # Tablas de resultados exportadas
│
└── README.md
```

---

## Metodología

### 1. Gestión de microdatos
Carga, limpieza y estructuración de los ficheros de microdatos del INE. Tratamiento de valores ausentes, codificación de variables y construcción del dataset analítico final.

### 2. Análisis Exploratorio (EDA)
Estadística descriptiva completa por grupos de gasto (medias, desviaciones, distribuciones). Visualización de distribuciones y correlaciones entre categorías de consumo.

### 3. Tests de hipótesis multivariantes
Validación estadística de diferencias significativas entre grupos mediante tests multivariantes (MANOVA, test de Hotelling). Comprobación de supuestos (normalidad, homocedasticidad).

### 4. Reducción de dimensionalidad
- **PCA** (Análisis de Componentes Principales): identificación de las dimensiones que explican la mayor varianza en el gasto de los hogares.
- **MDS** (Escalado Multidimensional): representación de las similitudes entre hogares en un espacio de baja dimensión.

### 5. Segmentación de consumidores
- **k-means clustering**: partición óptima de hogares con selección del número de clusters mediante criterio del codo y silueta.
- **Clustering jerárquico**: dendrograma para validar la estructura de grupos encontrada.
- Perfilado de cada segmento con interpretación orientada a decisiones de negocio.

---

## Resultados principales

- Se identificaron **X clusters** con perfiles de consumo claramente diferenciados (alta/media/baja intensidad de gasto en ocio, alimentación y vivienda).
- Las primeras **2 componentes principales** explican aproximadamente el **XX% de la varianza total**.
- Se encontraron diferencias estadísticamente significativas (p < 0.05) entre grupos en las principales categorías de gasto.

> *Nota: los resultados exactos están disponibles en la carpeta `outputs/`.*

---

## Cómo reproducir el análisis

### Requisitos

```r
# Instalar paquetes necesarios
install.packages(c(
  "tidyverse",   # manipulación y visualización de datos
  "FactoMineR",  # PCA y MDS
  "factoextra",  # visualización de análisis factorial
  "cluster",     # algoritmos de clustering
  "ggplot2",     # gráficos
  "corrplot"     # matrices de correlación
))
```

### Datos

Los microdatos de la EPF 2024 están disponibles de forma pública en el sitio web del INE:  
[https://www.ine.es/dyngs/INEbase/es/operacion.htm?c=Estadistica_C&cid=1254736176806](https://www.ine.es/dyngs/INEbase/es/operacion.htm?c=Estadistica_C&cid=1254736176806)

Descarga los ficheros y colócalos en `data/raw/` antes de ejecutar los scripts.

### Ejecución

Ejecuta los scripts en orden numérico desde RStudio o desde la consola de R:

```r
source("R/01_limpieza.R")
source("R/02_eda.R")
source("R/03_hipotesis.R")
source("R/04_pca_mds.R")
source("R/05_clustering.R")
```

---

## Autor

**An Fu Pham He**  
Estudiante de Estadística y Empresa · Universidad Carlos III de Madrid  
[LinkedIn](https://www.linkedin.com/in/an-fu-pham-he-4335732b8/) · anfuph07@gmail.com

---

## Licencia

Proyecto académico con fines educativos. Los microdatos utilizados son propiedad del INE y están sujetos a su política de uso.

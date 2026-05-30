# Análisis Multivariante del Consumo de los Hogares Españoles (EPF 2024)

## 📌 Descripción del Proyecto
Este repositorio contiene un análisis estadístico multivariante integral basado en los microdatos de la **Encuesta de Presupuestos Familiares (EPF) 2024** del Instituto Nacional de Estadística (INE). El proyecto aborda desde el preprocesamiento técnico de datos económicos hasta la creación de perfiles socioeconómicos complejos utilizando técnicas avanzadas de reducción de dimensionalidad y segmentación.

## 📊 Conjunto de Datos
Se ha trabajado con una muestra de **1000 hogares**, utilizando la clasificación oficial COICOP 2018 para las partidas de gasto. El dataset integra variables de naturaleza mixta:
* **6 Cuantitativas (Gasto en miles de €):** Comida, Vivienda, Sanidad, Ocio, Restaurantes/Hoteles y Gasto Total. 
* **2 Binarias:** Condición de asalariado y Otros ingresos.
* **3 Categóricas:** Nivel de estudios del sustentador, Tamaño del municipio y Tamaño del hogar.

## 🛠️ Tecnologías y Herramientas
* **Lenguaje:** R / RStudio.
* **Técnicas Estadísticas:** T² de Hotelling, MANOVA (Lambda de Wilks), PCA, MDS (Distancia G-Gower) y Clasificación Jerárquica.
* **Librerías principales:** `GGally`, `corrplot`, `cluster`, `ggplot2` y funciones de álgebra lineal para descomposición espectral.

---

## 🚀 Metodología y Fases del Análisis

### 1. Preprocesamiento y Contrastes de Hipótesis
* **Transformación Logarítmica:** Para corregir la asimetría positiva y la heterocedasticidad propias de los datos de gasto, se aplicó la transformación $log(x+1)$, permitiendo trabajar bajo supuestos de normalidad multivariante.
* **Validación Estadística:** Mediante los tests de **T² de Hotelling** y **Lambda de Wilks**, se confirmó que la estructura de gasto no es uniforme, sino que varía significativamente según la situación laboral y el nivel educativo ($p < 0.001$).

### 2. Análisis de Componentes Principales (PCA)
Se redujo la dimensionalidad de las variables de gasto, identificando dos dimensiones que explican el **64.76% de la varianza**:
* **Dimensión de Volumen (CP1):** Representa la capacidad de gasto total o nivel de riqueza del hogar.
* **Dimensión de Estilo de Vida (CP2):** Contrapone el gasto en necesidades básicas (Vivienda/Alimentación) frente al gasto discrecional (Ocio/Hostelería).

### 3. Escalamiento Multidimensional (MDS) y Datos Mixtos
Dada la presencia de variables cualitativas y cuantitativas, se empleó la **métrica G-Gower**. El análisis MDS permitió proyectar la proximidad entre hogares en un espacio bidimensional, revelando que el factor sociolaboral (ser asalariado o no) es el principal eje de diferenciación en la muestra.

### 4. Segmentación y Perfiles Socioeconómicos
Mediante clasificación jerárquica, se identificaron **4 arquetipos de consumo** en España:
1. **Consumo Discrecional (Clase Media-Alta):** Perfiles urbanos con alta inversión en ocio y servicios.
2. **Perfil de Vulnerabilidad:** Hogares con gasto bloqueado en necesidades de supervivencia y baja elasticidad.
3. **Estabilidad Familiar:** Hogares tradicionales de tamaño medio con consumo equilibrado.
4. **"Generación Inquilina":** Jóvenes altamente cualificados con alto gasto en servicios urbanos pero dificultad de ahorro debido a costes residenciales.

---

## 💡 Conclusiones de Negocio
El análisis revela que el consumo en España está fracturado estructuralmente por el ciclo de vida familiar y el nivel educativo. Se observa que el gasto en hostelería y cultura actúa como un fuerte indicador de estatus, mientras que el coste de la vivienda condiciona el resto de las decisiones financieras, especialmente en los perfiles más jóvenes y cualificados.

## 📂 Estructura del Repositorio
* `codigo_analisis/`: Script dividido en 3 partes con la limpieza, contrastes y algoritmos de reducción.
* `graficos/`: Carpetas de gráficos resultantes del análisis realizado.
* `memoria_tecnica/`: Documentación detallada de los resultados y la teoría aplicada. Estructurado en 3 partes.

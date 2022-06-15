
#* Este código ha sido creado y publicado originalmente en el idioma inglés,
#* como parte de un artículo científico titulado: 
#* 'An Introductory Tutorial to Cohort State-Transition Models in R for 
#* Cost-Effectiveness Analysis' 
#* Los autores originales son: 
#* - Fernando Alarid-Escudero <fernando.alarid@cide.edu>
#* - Eline Krijkamp
#* - Eva A. Enns
#* - Alan Yang
#* - M.G. Myriam Hunink
#* - Petros Pechlivanoglou
#* - Hawre Jalal
#* Por favor, cite el artículo si utiliza éste código.


#Traducción al Español: 
# Federico Rodriguez Cairoli, MD (IECS - R-HTA,LMICS)


#* A continuación se detallan los comandos necesarios para correr un modelo
#* de Markov hipotético (modelo de tipo cohorte con transición de estados) 
#* con variables que no varían a través del tiempo (tiempo independientes).
#* 
#* Ejemplo del caso: "Sick-sicker model".
#* 
#* "Sick-sicker model":
#* Este modelo se compone de cuatro estados, definidos por la condición de 
#* salud de una cohorte de pacientes hipotéticos: estado "Sano",
#* estado "Enfermo", estado "Enfermo-Severo", y estado "Muerte".
#* 
#* Las estrategias que se modelarán son: 
##* 1- Estándar de cuidado (en inglés: Standard of care, Soc), que representa 
#* la mejor atención disponible al momento, para los pacientes con esta 
#* enfermedad.
#* 
##* 2- Estrategia A (en inglés: Strategy A), que representa a la Terapia A, la 
#* cual es dada a todos los pacientes que se encuentren tanto en el estadio 
#* "Enfermo" (Sick) como en el estadio "Enfermo Severo" (Sicker). 
#* Sin embargo, esta terapia sólo mejora la calidad de vida de los pacientes 
#* en el estado "Enfermo" (Sick).
#* 
##* 3- Estragegia B (Strategy B), que representa a la Terapia B, la cual es dada
#* a todos los pacientes que se encuentren tantoen estadio "Enfermo" (Sick), 
#*  como en el estadio "Enfermo Severo" (Sicker). Esta terapia reduce la 
#*  progresión desde el estadio "Enfermo" (Sick) al estadio "Enfermo Severo" 
#*  (Sicker).
#*  
##* 4- Estrategia AB (Strategy AB), la cual representa una combinación de ambas 
#* terapias (A y B). Dada esta combinación, la progresión de la enfermedad desde
#*  el estadio "Enfermo" (Sick) al estadio "Enfermo Severo" (Sicker).es 
#* reducida y hay una mejora en la calidad de vida de los pacientes en el 
#* estado "Enfermo" (Sick).

#### 1. CONFIGURACION GENERAL ---- 
rm(list = ls())   # útil para remover cualquier variable en la memoria de R, 
#* y correr el modelo con mayor velocidad.


### Instalar los siguientes paquetes (sólo debe realizar esto UNA UNICA VEZ) 
install.packages("dplyr")      # importante para manejo de datos 
install.packages("tidyr")      # importante para manejo de datos 
install.packages("reshape2")   # importante para manejo de datos 
install.packages("ggplot2")    # importante para visualizar datos 
install.packages("ggrepel")    # importante para visualizar datos 
install.packages("ellipse")    # importante para visualizar datos 
install.packages("scales")     # Para signo dolar y escalas numéricas
install.packages("dampack")    # Para cálculo de ICER, entre otras funciones
install.packages("devtools")   # Para instalar paquetes desde GitHub
install.packages("rlang")      # Instala o actualiza este paquete
devtools::install_github("DARTH-git/darthtools") # para instalar paquete 
#darthtools desde GitHub 
install.packages("doParallel") # Para optimizar la corrida 
                              # de cualquier análisis de sensibilidad 

## Cargar los paquetes (Correr estos comandos cada vez que cierra y abre R) 
library(dplyr)
library(tidyr)
library(reshape2)  
library(ggplot2)    
library(ggrepel)   
library(ellipse)    
library(scales)     
library(dampack)  
library(darthtools)
library(doParallel) 

# Activar funciones necesarias, incluidas en carpeta "R" de este proyecto
source("R/Functions.R")

#### 2. PARAMETROS DEL MODELO ----

### Parámetros generales ----

cycle_length <- 1       # Duración del ciclo del modelo en años (1 año)
#(usar 1/12 en caso de ser ciclo mensual)
n_age_init <- 25        # Edad de los pacientes al inicio (25 años)
n_age_max  <- 100       # Máxima edad de seguimiento de los pacientes(100 años)
n_cycles <- (n_age_max - n_age_init)/cycle_length # Número total de ciclos 
# en el modelo (Time horizon)
v_names_states <- c("H",  # Los cuatro estados del modelo:
                    "S1", # Sano (Healthy: "H"), Enfermo (Sick:"S1"), 
                    # Enfermo Severo (Sicker: "S2"), Muerto (Dead:"D")
                    "S2", 
                    "D") 

n_states <- length(v_names_states)     # Número total de estados

### Tasa de descuento anual para costos y QALYs ----
d_c <- 0.03 # Tasa de descuento anual para costos (3%) 
d_e <- 0.03 # Tasa de descuento anual para QALYs (3%)  

### Estrategias ----
v_names_str <- c("Standard of care",    # Guardar el nombre de las estrategias
                 "Strategy A", 
                 "Strategy B",
                 "Strategy AB") 
n_str       <- length(v_names_str)        # Número total de las estrategias

### Corrección de ciclo ----
v_wcc <- gen_wcc(n_cycles = n_cycles,  # Corrección de ciclo anual  usando regla 
                                       # Simpson 1/3 
                 method = "Simpson1/3") 

### Tasas anuales de transición y Hazard ratios (HRs) ----
r_HD   <- 0.002 # Tasa anual de mortalidad por todas las causas (constante) en 
                #  sujetos en estado "Sano".
r_HS1  <- 0.15  # Tasa anual de volverse "Enfermo" en sujetos en estado "Sano" 
r_S1H  <- 0.5   # Tasa anual de volverse "Sano" en sujetos en estado "Enfermo" 
r_S1S2 <- 0.105 # Tasa anual de volverse "Enfermo severo" en sujetos en 
                # estado "Enfermo" 
hr_S1  <- 3     # Hazard ratio de muerte en "Enfermos" en comparación a "Sanos" 
hr_S2  <- 10    # Hazard ratio de muerte en "Enfermos severos" en comparación 
                # a "Sanos" 

### Efectividad del tratamiento B ----
hr_S1S2_trtB <- 0.6  # Hazard ratio de volverse "Enfermo severo" en sujetos en 
                     # estado "Enfermo" bajo el tratamiento B

#### Costos y utilidades 
### Costos ----
c_H    <- 2000  # Costo anual por ser "Sano (H)"
c_S1   <- 4000  # Costo anual por estar "Enfermo (S1)"
c_S2   <- 15000 # Costo anual por ser "Enfermo Severo (S2)"
c_D    <- 0     # Costo anual en estado de "Muerte (D)"
c_trtA <- 12000 # Costo anual por recbibir Terapia A
c_trtB <- 13000 # Costo anual por recbibir Terapia B
### Utilidades ----
u_H    <- 1     # Utilidad anual por ser "Sano (H)"
u_S1   <- 0.75  # Utilidad anual por estar "Enfermo (S1)"
u_S2   <- 0.5   # Utilidad anual por ser "Enfermo Severo (S2)"
u_D    <- 0     # Utilidad anual en estado de "Muerte (D)"

### Efectividad del tratamiento A ----
u_trtA <- 0.95  # Utilidad anual en caso de recbibir Terapia A

### Descuento ponderado (en caso de querer aplicarlo) para costos y QALYs 
v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))

#### 3. PROCESAMIENTO DE LOS PARAMETROS ----

### Aplicar Hazard ratios ----
##* Calcular las tasas anuales de mortalidad para los estados "Enfermo" 
#* y "Enfermo Severo"
r_S1D <- r_HD * hr_S1 # Tasa anual de mortalidad para "Enfermos (S1)"
r_S2D <- r_HD * hr_S2 # Tasa anual de mortalidad para "Enfermos Severos (S2)"

### Probabilidades de Transición ----
#* Convertir tasas anuales a probabilidades de transición
#* Función "rate_to_prob" del paquete `darthtools` 
p_HS1  <- rate_to_prob(r = r_HS1, t = cycle_length) # Probabilidad anual de 
# volverse "Enfermo (S1)" estando en estadio "Sano (H)" 
p_S1H  <- rate_to_prob(r = r_S1H, t = cycle_length) # Probabilidad anual de 
# volverse "Sano (H)" estando en estadio "Enfermo (S1)"
p_S1S2 <- rate_to_prob(r = r_S1S2, t = cycle_length)# Probabilidad anual de 
# volverse "Enfermo Severo (S2)" estando en estadio "Enfermo (S1)" 
p_HD   <- rate_to_prob(r = r_HD, t = cycle_length)  # Probabilidad anual de 
# morir estando en estadio "Sano (H)" 
p_S1D  <- rate_to_prob(r = r_S1D, t = cycle_length) # Probabilidad anual de 
# morir estando en estadio "Enfermo (S1)" 
p_S2D  <- rate_to_prob(r = r_S2D, t = cycle_length)  # Probabilidad anual de 
# morir estando en estadio "Enfermo Severo (S2)" 


### Eficacia de Terapia B ----
##Probabilidad de transición de pasar a "Enfermo Severo (S2)", 
# estando en estadio Enfermo (S1), bajo la Terapia B:
##* Primero, aplicar Hazar ratio a la tasa anual de pasar a estadio "Enfermo 
#* Severo (S2)", estando en estadio "Enfermo (S1), bajo la Terapia B
r_S1S2_trtB <- r_S1S2 * hr_S1S2_trtB
##* Luego transformar dicha tasa anual a probabilidad. 
p_S1S2_trtB <- rate_to_prob(r = r_S1S2_trtB, t = cycle_length) 

#### 4. GENERACION DE COHORTES Y MATRICES DE PROBABILIDADES DE TRANSICION ----
### Cohorte de pacientes - Vector inicial ----
#* Todos los pacientes comienzan estando "Sanos (H)". A modo de simplicidad, 
#* se utiliza un valor "1", que representa el 100% de la cohorte de pacientes
#*  al inicio (pudiera ser cualquier número de pacientes iniciales: 
#*  por ejemplo, 1, 100, 10000, etc. Sea 1, 100 o 10000, esto no tiene ninguna 
#*  repercusión en los resultados finales del ICER.
v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) 
v_m_init

### Cohortes para cada estrategias de análisis  ----
## Cohorte para "Estándard de cuidado (SOC)" ----
m_M <- matrix(NA, 
              nrow = (n_cycles +1), ncol = n_states, 
              dimnames = list(0:n_cycles, v_names_states))
#* Anclar el vector inicial como primera fila de nuestra cohorte (matriz)
m_M[1, ] <- v_m_init

## Cohorte para Estrategias A, B, y AB ----
#* La estructura inicial es igual a la creada para "Estándard de cuidado (SOC)"
m_M_strA  <- m_M # Estrategia A (Terapia A)
m_M_strB  <- m_M # Estrategia B (Terapia B)
m_M_strAB <- m_M # Estrategia AB (Terapia C)

### Matrices de probabilidades de transición ---- 
## Matriz inicial para "Estándard de cuidado (SOC)"   
 m_P <- matrix(0, 
              nrow = n_states, ncol = n_states, 
              dimnames = list(v_names_states, 
                              v_names_states)) 
m_P

#* Completar la matriz con las probabilidades de transición correspondientes
## Matriz de probabilidad de transición para "Estándard de cuidado (SOC)" ----
#* Desde "Sano (H)"
m_P["H", "H"]   <- (1 - p_HD) * (1 - p_HS1)
m_P["H", "S1"]  <- (1 - p_HD) * p_HS1 
m_P["H", "D"]   <- p_HD
#* Desde estadio "Enfermo (S1)"
m_P["S1", "H"]  <- (1 - p_S1D) * p_S1H
m_P["S1", "S1"] <- (1 - p_S1D) * (1 - (p_S1H + p_S1S2))
m_P["S1", "S2"] <- (1 - p_S1D) * p_S1S2
m_P["S1", "D"]  <- p_S1D
#* Desde estadio "Enfermo Severo (S2)"
m_P["S2", "S2"] <- 1 - p_S2D
m_P["S2", "D"]  <- p_S2D
#* Desde estadio "Muerte (D)"
m_P["D", "D"]   <- 1

## Matriz de probabilidades de transición para "Estrategia A (terapia A)"   -- 
# Dado que las probabilidades son las mismas, crear una copia de matriz  
# creada para "Estándard de cuidado  (SoC)" ----
m_P_strA <- m_P

## Matriz de probabilidades de transición para "Estrategia B (terapia B)" -- 
m_P_strB <- m_P
#* Modificar SOLO las probabilidades de transición que involucran estadios 
#* de "Enfermos (S1)" y "Enfermos Severos (S2) 
m_P_strB["S1", "S1"] <- (1 - p_S1D) * (1 - (p_S1H + p_S1S2_trtB))
m_P_strB["S1", "S2"] <- (1 - p_S1D) * p_S1S2_trtB

## Matriz de probabilidades de transición para "Estrategia AB (terapias A y B)"  
#* Crear una copia de Matriz para Estrategia B, dado que es la misma matriz 
m_P_strAB <- m_P_strB

## Validación de probabildiades de transición [0, 1] ----
check_transition_probability(m_P,      verbose = TRUE)  # m_P >= 0 && m_P <= 1
check_transition_probability(m_P_strA, verbose = TRUE)  # m_P_strA >= 0 && m_P_strA <= 1
check_transition_probability(m_P_strB, verbose = TRUE)  # m_P_strB >= 0 && m_P_strB <= 1
check_transition_probability(m_P_strAB, verbose = TRUE) # m_P_strAB >= 0 && m_P_strAB <= 1
### Chequear que todas las filas sumen 1 ----
check_sum_of_transition_array(m_P,      n_states = n_states, n_cycles = n_cycles, verbose = TRUE)  # Suma de filas (m_P) == 1
check_sum_of_transition_array(m_P_strA, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)  # Suma de filas (m_P_strA) == 1
check_sum_of_transition_array(m_P_strB, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)  # Suma de filas (m_P_strB) == 1
check_sum_of_transition_array(m_P_strAB, n_states = n_states, n_cycles = n_cycles, verbose = TRUE) # Suma de filas (m_P_strAB) == 1


####  5. CORRER EL MODELO  ----
## Integrar matriz de probabilidades de transición con matriz de cohorte ----
for(t in 1:n_cycles){
  # Para "Estandar of care (SoC)"
  m_M[t + 1, ] <- m_M[t, ] %*% m_P
  # Para "Estrategia A (Terapia A)"
  m_M_strA[t + 1, ] <- m_M_strA[t, ] %*% m_P_strA
  # Para "Estrategia B (Terapia B)"
  m_M_strB[t + 1, ] <- m_M_strB[t, ] %*% m_P_strB
  # Para "Estrategia AB (Terapia A y B)"
  m_M_strAB[t + 1, ] <- m_M_strAB[t, ] %*% m_P_strAB
}

## Visualizar las cohortes definitivas
m_M
m_M_strA
m_M_strB
m_M_strAB


## Guardar las cohortes completas en una lista ----
l_m_M <- list(m_M,
              m_M_strA,
              m_M_strB,
              m_M_strAB)
names(l_m_M) <- v_names_str

### Visualización de las cohortes entre estados en el tiempo de análisis ----
##* Gráfico EJEMPLO para "Estandar de cuidado (SoC)"
plot_trace( m_M) 

### Integrar Utilidades y costos ----
##* Vector de utilidades por estado para "Estandar of care (SoC)"
v_u_SoC    <- c(H  = u_H, 
                S1 = u_S1, 
                S2 = u_S2, 
                D  = u_D) * cycle_length
##* Vector Vector de costos por estado para "Estandar of care (SoC)"
v_c_SoC    <- c(H  = c_H, 
                S1 = c_S1,
                S2 = c_S2, 
                D  = c_D) * cycle_length
##* Vector de utilidades por estado para "Estrategia A"
v_u_strA   <- c(H  = u_H, 
                S1 = u_trtA, 
                S2 = u_S2, 
                D  = u_D) * cycle_length
##* Vector de costos por estado para "Estrategia A"
v_c_strA   <- c(H  = c_H, 
                S1 = c_S1 + c_trtA,
                S2 = c_S2 + c_trtA, 
                D  = c_D)
##* Vector de utilidades por estado para "Estrategia B"
v_u_strB   <- c(H  = u_H, 
                S1 = u_S1, 
                S2 = u_S2, 
                D  = u_D) * cycle_length
##* Vector de costos por estado para "Estrategia B"
v_c_strB   <- c(H  = c_H, 
                S1 = c_S1 + c_trtB, 
                S2 = c_S2 + c_trtB, 
                D  = c_D) * cycle_length
##* Vector de utilidades por estado para "Estrategia AB"
v_u_strAB  <- c(H  = u_H, 
                S1 = u_trtA, 
                S2 = u_S2, 
                D  = u_D) * cycle_length
##* Vector de costos por estado para "Estrategia AB"
v_c_strAB  <- c(H  = c_H, 
                S1 = c_S1 + (c_trtA + c_trtB), 
                S2 = c_S2 + (c_trtA + c_trtB), 
                D  = c_D) * cycle_length


##* Almacenar en una lista los vectores de las utilidades para todas 
#* las estrategias de análisis
l_u   <- list(SQ = v_u_SoC,
              A  = v_u_strA,
              B  = v_u_strB,
              AB = v_u_strAB)
##* Almacenar en una lista los vectores de los costos para todas 
#* las estrategias de análisis
l_c   <- list(SQ = v_c_SoC,
              A  = v_c_strA,
              B  = v_c_strB,
              AB = v_c_strAB)

#* Asignar los nombres completos de cada estrategia a las listas recién creadas
names(l_u) <- names(l_c) <- v_names_str

#### 6. RESULTADOS ----
##* Para obtener resultados, primero crear vectores vacíos con el objetivo
##* de almacenar todas las utilidades y costos totales por estrategia  
v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str

## Segundo, crear loop para calcular utilidades y costos totales  
for (i in 1:n_str) {
  v_u_str <- l_u[[i]]   
  v_c_str <- l_c[[i]]  
  
##* Qalys y costos esperados por ciclo 
  v_qaly_str <- l_m_M[[i]] %*% v_u_str # Suma de utilidades de todos los estados
                                       # por ciclo
  v_cost_str <- l_m_M[[i]] %*% v_c_str # Suma de costos de todos los estados
                                       # por ciclo
  ###* QALYs y costos totales descontados por estrategia (con ajuste por ciclo)
  #* QALYs
  v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
  #* Costos
  v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
}

# Análisis de costo-efectividad (CEA) ----
## Razón de costo-efectividad incremental (ICERs) ----
#* La última versión de la función que se utilizará a continuación 
#* "calculate_icers", se encuentra disponible en el paquete "dampack"
df_cea <- dampack:: calculate_icers(cost       = v_tot_cost, 
                          effect     = v_tot_qaly,
                          strategies= v_names_str)
df_cea

## Tabla CEA óptima  ----
table_cea <- format_table_cea(df_cea) 

## Frontera de costo-efectividad -----
#* La última versión de la función que se utilizará a continuación 
#* se encuentra disponible en el paquete "dampack".
plot(df_cea, label = "all", txtsize = 16) +
   expand_limits(x = max(table_cea$QALYs) + 0.1) +
   theme(legend.position = c(0.8, 0.2))





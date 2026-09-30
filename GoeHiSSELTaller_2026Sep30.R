#########################################################################
#########################################################################
# Taller GeoHiSEE

#########################################################################
#########################################################################

#########################################################################
# 1) Instalación y carga de paquetes
install.packages("hisse")
# Instalar devtools si no lo tienes
install.packages("devtools")
# Instalar hisse desde GitHub
library(devtools)
install_github(repo = "thej022214/hisse", ref = "master")
# remotes::install_github("thej022214/hisse")
# pak::pak("thej022214/hisse")

# Cargar los paquetes necesarios
suppressWarnings(library(hisse))
suppressWarnings(library(diversitree))
library(ape)
library(phytools)

library(readxl)

##########################################################################
# 2. Cargar y verificar el árbol en R

# Cargar el árbol de TimeTree

phy <- read.tree("inga_species.nwk") # también puede ser .tree
# Verificar propiedades
print(phy) # resumen general
phy$tip.label # deben aparecer nombres de las especies
is.rooted(phy) # para ver si es enraizado, debe ser TRUE
is.binary(phy) # preferible si está completamente resuelto, FALSE es normal en
# árboles unsmoothed
is.ultrametric(phy) # Debe serlo, es decir la longitud de las ramas ajustadas al
# tiempo absoluto, se espera TRUE


# Verificar magnitud del problema si no es ultramétrico
distancias <- dist.nodes(phy)[1:length(phy$tip.label), 1]
range(distancias) # si el rango es < 1e-6, es solo precisión numérica
# Corregir extendiendo las ramas terminales
# phy <- force.ultrametric(phy, method = "extend")
is.ultrametric(phy) # ahora TRUE

###########################################################################
# 3. Preparar los datos de rangos geográficos


# rangos geográficos para las especies de Inga en la Mata Atlántica y Amazonas

# leer tabla de la distribución del género Inga
location.inga <- read_excel("Location_Inga.xlsx", sheet = "Table008 (Page 6-17)")
head(location.inga)
View(location.inga)

#Resumen de la distribución de las especies del género Inga en las dos regiones a comparar (Adaptado de Nicholls,
#James A., Jens J. Ringelberg, Kyle G. Dexter, et al. «Continuous colonization of the Atlantic coastal rain
#forests of South America from Amazônia». Proceedings of the Royal Society B: Biological Sciences 292, n.º 2039
#(2025): 20241559. https://doi.org/10.1098/rspb.2024.1559.)

regions <- table(location.inga$species, location.inga$region)
rownames(regions)

#Ajustar los nombre de las especies taxonómicas para que coincida a los del árbol filogenético
rownames(regions) <- sub("\r\n", "_", rownames(regions))

#Verificar si coinciden
phy$tip.label %in% rownames(regions)

#Las que coinciden
phy$tip.label[which(phy$tip.label %in% rownames(regions))]
length(phy$tip.label[which(phy$tip.label %in% rownames(regions))])
# Las que no
phy$tip.label[which(!phy$tip.label %in% rownames(regions))]

# filtrar solo las que sí coinciden, por ahora
regions.df <- as.data.frame.matrix(regions)
dim(regions.df)
rownames(regions.df)
regions.phy <- regions.df[rownames(regions.df) %in% phy$tip.label, c(1, 3)]
colnames(regions.phy)
dim(regions.phy)
regions.phy$range <- rep(NA, nrow(regions.phy))
# si endémico de Amazonas
regions.phy$range[which(regions.phy$`Mata Atlantica` == 0 &
                                regions.phy$Amazon != 0)] <- 1
# si es endémico de mata Atlántica
regions.phy$range[which(regions.phy$`Mata Atlantica` != 0 &
                                regions.phy$Amazon == 0)] <- 2

# si es de amplia distribución
regions.phy$range[which(regions.phy$`Mata Atlantica` != 0 &
                                regions.phy$Amazon != 0)] <- 0
regions.phy[is.na(regions.phy[, 3]), ]

#agregar algunas distribuciones que no se encuentran en Nicholls, James A., Jens J. Ringelberg, Kyle G. Dexter,
#et al. «Continuous colonization of the Atlantic coastal rain forests of South America from Amazônia».
#Proceedings of the Royal Society B: Biological Sciences 292, n.º 2039 (2025): 20241559.
#https://doi.org/10.1098/rspb.2024.1559.

spp <- sort(phy$tip.label[which(!phy$tip.label %in% rownames(regions))])
ranges <- c(0, 0, 0, 1, 1, 1, 0, rep(NA, 23))
length(spp)
length(ranges)
regions.phy.df <- data.frame(especie = c(rownames(regions.phy), spp),
                             rangos = c(regions.phy$range, ranges))
dim(regions.phy.df)
length(phy$tip.label)

unique(regions.phy.df$especie)
#regions.phy.df coinciden los rangos con los tips del árbol filogenético?!phy$tip.label %in% regions.phy.df$especie
length(which(!phy$tip.label %in% regions.phy.df$especie))


# Verificar coincidencias
matched <- regions.phy.df$especie %in% phy$tip.label
cat("Especies coincidentes:",
    sum(matched),
    "de",
    nrow(regions.phy.df),
    "\n")


# Puesto que las filogenias usualmente no son completas, aqui es necesario “decirle” al
# programa cual es la proporción de especies muetreadas en cada una de las
# distribuciones. Para esto definimos el vector f.
# Vamos a suponer que el muestreo es completo, pero lo pueden modificar a partir de la
# consulta del muestreo taxonómico en su árbol.

f <- c(1, 1, 1)


# como hay especies que no tienen información de la distribución las eliminaré de la filogenia
str(phy)
? drop.tip

phy.drop <- drop.tip(phy = phy, tip = regions.phy.df$especie[which(is.na(regions.phy.df$rangos))])

phy.drop$tip.label
length(phy.drop$tip.label)
length(regions.phy.df$especie[which(!is.na(regions.phy.df$rangos))])

regions.phy.df <- regions.phy.df[which(!is.na(regions.phy.df$rangos)), ]
# Verificar coincidencias
matched <- regions.phy.df$especie %in% phy.drop$tip.label
cat("Especies coincidentes:",
    sum(matched),
    "de",
    nrow(regions.phy.df),
    "\n")

#################################################################################
# 4. Construir la matriz de transición


# La matriz de transición representa las tasas de los posibles cambios de distribución. No
# indica cuántas especies hay en cada región, sino qué tan rápidamente el modelo
# permite que los linajes cambien de estado geográfico. De esta manera, la matriz define
# qué transiciones están permitidas y qué parámetros se comparten entre estados.

# Modelo GeoSSE sin estados ocultos
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits = 0)
# Modelo # Modelo # Modelo GeoHiSSE con 1 estado oculto
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits = 1)
# Modelo nulo (transiciones simplificadas)
trans.rate.null <- TransMatMakerGeoHiSSE(hidden.traits = 1, make.null = TRUE)

###############################################################################
# 5. Definir los vectores de tasas para los tres modelos

# Modelo 1: sin estados ocultos
turnover_geo <- c(1, 2, 3) # s0, s1, s01 libres
eps_geo <- c(1, 2) # x0, x1 libres
hidden_geo <- FALSE
trans_geo <- TransMatMakerGeoHiSSE(hidden.traits = 0)
# Modelo 2: con 1 estado oculto, todas las tasas libres
turnover_hid <- c(1, 2, 3, 4, 5, 6) # cada estado con su propia tasa
eps_hid <- c(1, 2, 3, 4) # cada estado con su propia fracción
hidden_hid <- TRUE
trans_hid <- TransMatMakerGeoHiSSE(hidden.traits = 1)

# Modelo 3: nulo (hidden states con tasas idénticas)
turnover_null <- c(1, 2, 3, 1, 2, 3) # s0A = s0B, s1A = s1B, s01A = s01B
eps_null <- c(1, 2, 1, 2) # x0A = x0B, x1A = x1B
hidden_null <- TRUE # sigue teniendo hidden states
trans_null <- TransMatMakerGeoHiSSE(hidden.traits = 1, make.null = TRUE)

##############################################################################
# 6. Ejecutar los tres modelos GeoHiSSE

## 6.1. GeoSSE estándar (sin estados ocultos)GeoSSE estándar (sin estados ocultos)

# Hipótesis: las tasas de diversificación dependen del rango geográfico, pero no existe
# heterogeneidad oculta.
mod1 <- GeoHiSSE(
        phy = phy.drop,
        data = regions.phy.df,
        f = f,
        turnover = turnover_geo,
        # c(1, 2, 3)
        eps = eps_geo,
        # c(1, 2)
        hidden.states = hidden_geo,
        # FALSE
        trans.rate = trans_geo,
        # TransMatMakerGeoHiSSE(hidden.traits = 0)
        sann = TRUE,
        sann.its = 1000
)
# Fit
# lnL              AIC             AICc           n.taxa n.hidden.classes
# -107.4543         228.9085         230.4220          82.0000           1.0000
#
# Model parameters:
#
#         tau00A    tau11A    tau01A     ef00A     ef11A  d00A_01A  d11A_01A
# 1.166552 17.083608  7.019174  3.000000  1.539356  1.064936 16.779788


## 6.2. GeoHiSSE con un estado oculto y tasas libres

# Hipótesis: además de la dependencia del rango geográfico, existe heterogeneidad
# oculta en las tasas que no se explica por la geografía.

mod2 <- GeoHiSSE(
        phy = phy.drop,
        data = regions.phy.df,
        f = f,
        turnover = turnover_hid,
        # c(1, 2, 3, 4, 5, 6)
        eps = eps_hid,
        # c(1, 2, 3, 4)
        hidden.states = hidden_hid,
        # TRUE
        trans.rate = trans_hid,
        # TransMatMakerGeoHiSSE(hidden.traits = 1)
        sann = TRUE,
        sann.its = 1000
)

# Fit
# lnL              AIC             AICc           n.taxa n.hidden.classes
# -105.8874         241.7747         249.0474          82.0000           2.0000
#
# Model parameters:
#
#         tau00A       tau11A       tau01A        ef00A        ef11A     d00A_01A
# 2.071279e-09 9.605345e-01 2.401338e-01 5.678671e-07 3.000000e+00 2.128658e-09
# d11A_01A     d00A_00B     d11A_11B     d01A_01B       tau00B       tau11B
# 1.000000e+02 1.345857e+00 1.345857e+00 1.345857e+00 3.701405e+00 2.141278e-09
# tau01B        ef00B        ef11B     d00B_01B     d11B_01B     d00B_00A
# 2.449432e+00 5.111279e-01 2.996020e+00 2.756927e-01 3.918563e-09 1.345857e+00
# d11B_11A     d01B_01A
# 1.345857e+00 1.345857e+00


## 6.3. Modelo nulo (estados ocultos sin efecto)

#Hipótesis: existen estados ocultos, pero no cambian las tasas. Es decir, la
# heterogeneidad oculta que el modelo 2 "captura" es indistinguible del ruido.

# Modelo nulo (estados ocultos sin efecto)
# Hipótesis: existen estados ocultos, pero no cambian las tasas. Es decir, la
# heterogeneidad oculta que el modelo 2 "captura" es indistinguible del ruido.

mod3 <- GeoHiSSE(
        phy = phy.drop,
        data = regions.phy.df,
        f = f,
        turnover = turnover_null,
        # c(1, 2, 3, 1, 2, 3)
        eps = eps_null,
        # c(1, 2, 1, 2)
        
        hidden.states = hidden_null,
        # TRUE (sigue teniendo estados ocultos)
        trans.rate = trans_null,
        # TransMatMakerGeoHiSSE(hidden.traits = 1, make.null
        sann = TRUE,
        sann.its = 1000
)
# Fit
# lnL              AIC             AICc           n.taxa n.hidden.classes
# -107.5023         231.0046         232.9772          82.0000           2.0000
#
# Model parameters:
#
#         tau00A       tau11A       tau01A        ef00A        ef11A     d00A_01A
# 1.130944e+00 1.635051e+01 6.512240e+00 3.000000e+00 1.624689e+00 1.100245e+00
# d11A_01A     d00A_00B     d11A_11B     d01A_01B       tau00B       tau11B
# 2.066411e+01 3.099127e-08 3.099127e-08 3.099127e-08 1.130944e+00 1.635051e+01
# tau01B        ef00B        ef11B     d00B_01B     d11B_01B     d00B_00A
# 6.512240e+00 3.000000e+00 1.624689e+00 1.100245e+00 2.066411e+01 3.099127e-08
# d11B_11A     d01B_01A
# 3.099127e-08 3.099127e-08

# Revisar problemas...
# Log-verosimilitudes (deben ser negativas)
mod1$loglik
mod2$loglik
mod3$loglik
# AIC y AICc
mod1$AIC
mod2$AIC
mod3$AIC
# Verificar que ningún parámetro tocó los límites del espacio de búsqueda
# (valores como 100, 1000 o 1e-8 suelen indicar problemas de convergencia)
mod1$solution
mod2$solution
mod3$solution

########################################################################################
# 7. Comparar los tres modelos con AIC / AICc

# Comparación directa
GetAICWeights(list(
        model1 = mod1,
        model2 = mod2,
        model3 = mod3
), criterion = "AIC")
? AIC
# model1      model2      model3
# 0.739519565 0.001188764 0.259291671

# AICc para cada modelo
mod1$AICc
mod2$AICc
mod3$AICc

# O calcularlo manualmente si el paquete no lo devuelve
AICc <- function(mod, n) {
        k <- length(mod$solution)
        mod$AIC + (2 * k * (k + 1)) / (n - k - 1)
}
n <- length(phy$tip.label)
AICc_val <- c(AICc(mod1, n), AICc(mod2, n), AICc(mod3, n))
AICc_val #-843.5359 -830.6697 -841.4398

# Reunir AICc en un vector
aicc <- c(mod1$AICc, mod2$AICc, mod3$AICc)
names(aicc) <- c("sin_ocultos", "1_oculto_libre", "nulo")
# ΔAICc
delta <- aicc - min(aicc)
# Pesos AIC
pesos <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))

# Tabla resumen
tabla <- data.frame(
        modelo = names(aicc),
        k = c(5, 10, 5),
        # parámetros libres de cada modelo
        AICc = round(aicc, 2),
        delta_AICc = round(delta, 2),
        peso_AICc = round(pesos, 3),
        loglik = round(c(mod1$loglik, mod2$loglik, mod3$loglik), 2)
)
print(tabla)
#                               modelo  k   AICc delta_AICc peso_AICc  loglik
# sin_ocultos       sin_ocultos  5 230.42       0.00     0.782 -107.45
# 1_oculto_libre 1_oculto_libre 10 249.05      18.63     0.000 -105.89
# nulo                     nulo  5 232.98       2.56     0.218 -107.50

# Guardar los resultados
# write.csv(tabla, "comparacion_modelos.csv", row.names = FALSE)
# saveRDS(list(mod1 = mod1, mod2 = mod2, mod3 = mod3), "tres_modelos.rds")

#########################################################################################
# 8. Interpretar las tasas de diversificación

# Como dos modelos tuvieron mayor soporte empírico, se examinan los dos:

## 8.1 Mejor modelo1
mejor_modelo_1 <- mod1
# Extraer parámetros
mejor_modelo_1$solution
class(mejor_modelo_1$solution)
names(mejor_modelo_1$solution)

# Conversión manual de turnover y eps a lambda y mu
# Si tienes los valores de turnover (tau) y extinction fraction (eps):


# tasas de especiación (λ)
which(startsWith(names(mejor_modelo_1$solution), "tau"))
tau <- mejor_modelo_1$solution[which(startsWith(names(mejor_modelo_1$solution), "tau"))]
tau
length(tau)
tau[-seq(3, 30, 3)]
which(startsWith(names(mejor_modelo_1$solution), "ef"))
eps <- mejor_modelo_1$solution[which(startsWith(names(mejor_modelo_1$solution), "ef"))]
eps
lambda <- tau[-seq(3, 30, 3)] / (1 + eps)
lambda
names(lambda) <- sub("tau", "", names(lambda))
# tau00A       tau11A       tau00B       tau11B       tau00C       tau11C
# 0.1077349151 0.0001439461 0.8419878876 5.2262792342 0.0000000000 0.0000000000
# tau00D       tau11D       tau00E       tau11E       tau00F       tau11F
# 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
# tau00G       tau11G       tau00H       tau11H       tau00I       tau11I
# 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
# tau00J       tau11J
# 0.0000000000 0.0000000000



lambda.matrix <- matrix(lambda,
                        nrow = 10,
                        ncol = 2,
                        byrow = T)
rownames(lambda.matrix) <- LETTERS[1:10]
colnames(lambda.matrix) <-  c("00", "11")

lambda.matrix
# 00       11
# A 0.2916381 6.727536
# B 0.0000000 0.000000
# C 0.0000000 0.000000
# D 0.0000000 0.000000
# E 0.0000000 0.000000
# F 0.0000000 0.000000
# G 0.0000000 0.000000
# H 0.0000000 0.000000
# I 0.0000000 0.000000
# J 0.0000000 0.000000

#tasas de extinción (μ)

mu <- tau[-seq(3, 30, 3)] * eps / (1 + eps)
# tau00A       tau11A       tau00B       tau11B       tau00C       tau11C
# 3.232047e-01 3.773917e-13 2.029284e+00 4.488215e+00 0.000000e+00 0.000000e+00
# tau00D       tau11D       tau00E       tau11E       tau00F       tau11F
# 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
# tau00G       tau11G       tau00H       tau11H       tau00I       tau11I
# 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
# tau00J       tau11J
# 0.000000e+00 0.000000e+00

mu.matrix <- matrix(mu,
                    nrow = 10,
                    ncol = 2,
                    byrow = T)
rownames(mu.matrix) <- LETTERS[1:10]
colnames(mu.matrix) <-  c("00", "11")

mu.matrix
# 00       11
# A 0.8749142 10.35607
# B 0.0000000  0.00000
# C 0.0000000  0.00000
# D 0.0000000  0.00000
# E 0.0000000  0.00000
# F 0.0000000  0.00000
# G 0.0000000  0.00000
# H 0.0000000  0.00000
# I 0.0000000  0.00000
# J 0.0000000  0.00000


## 8.3 Mejor modelo3
mejor_modelo_3 <- mod3
# Extraer parámetros
mejor_modelo_3$solution
class(mejor_modelo_3$solution)
names(mejor_modelo_3$solution)

# Conversión manual de turnover y eps a lambda y mu
# Si tienes los valores de turnover (tau) y extinction fraction (eps):


# tasas de especiación (λ)
which(startsWith(names(mejor_modelo_3$solution), "tau"))
tau3 <- mejor_modelo_3$solution[which(startsWith(names(mejor_modelo_3$solution), "tau"))]
tau3
length(tau3)
tau3[-seq(3, 30, 3)]
which(startsWith(names(mejor_modelo_3$solution), "ef"))
eps3 <- mejor_modelo_3$solution[which(startsWith(names(mejor_modelo_3$solution), "ef"))]
eps3
lambda3 <- tau3[-seq(3, 30, 3)] / (1 + eps3)
lambda3
names(lambda3) <- sub("tau", "", names(lambda3))
# tau00A   tau11A   tau00B   tau11B   tau00C   tau11C   tau00D   tau11D   tau00E
# 0.282736 6.229503 0.282736 6.229503 0.000000 0.000000 0.000000 0.000000 0.000000
# tau11E   tau00F   tau11F   tau00G   tau11G   tau00H   tau11H   tau00I   tau11I
# 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
# tau00J   tau11J
# 0.000000 0.000000

lambda.matrix3 <- matrix(lambda3,
                         nrow = 10,
                         ncol = 2,
                         byrow = T)
rownames(lambda.matrix3) <- LETTERS[1:10]
colnames(lambda.matrix3) <-  c("00", "11")

lambda.matrix3
# 00       11
# A 0.282736 6.229503
# B 0.282736 6.229503
# C 0.000000 0.000000
# D 0.000000 0.000000
# E 0.000000 0.000000
# F 0.000000 0.000000
# G 0.000000 0.000000
# H 0.000000 0.000000
# I 0.000000 0.000000
# J 0.000000 0.000000

#tasas de extinción (μ)

mu3 <- tau3[-seq(3, 30, 3)] * eps3 / (1 + eps3)
# tau00A     tau11A     tau00B     tau11B     tau00C     tau11C     tau00D
# 0.8482079 10.1210052  0.8482079 10.1210052  0.0000000  0.0000000  0.0000000
# tau11D     tau00E     tau11E     tau00F     tau11F     tau00G     tau11G
# 0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
# tau00H     tau11H     tau00I     tau11I     tau00J     tau11J
# 0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000

mu.matrix3 <- matrix(mu3,
                     nrow = 10,
                     ncol = 2,
                     byrow = T)
rownames(mu.matrix3) <- LETTERS[1:10]
colnames(mu.matrix3) <-  c("00", "11")

mu.matrix3
# 00       11
# A 0.8482079 10.12101
# B 0.8482079 10.12101
# C 0.0000000  0.00000
# D 0.0000000  0.00000
# E 0.0000000  0.00000
# F 0.0000000  0.00000
# G 0.0000000  0.00000
# H 0.0000000  0.00000
# I 0.0000000  0.00000
# J 0.0000000  0.00000

# tasas de dispersión
dispersion <- mejor_modelo_3$solution[which(startsWith(names(mejor_modelo_3$solution), "d"))]
length(dispersion)
sort(dispersion)
sort(dispersion[which(dispersion>3.099127e-08)], decreasing = T)
# d11A_01A     d11B_01B     d00A_01A     d00B_01B     d00A_00B     d11A_11B 
# 2.066411e+01 2.066411e+01 1.100245e+00 1.100245e+00 3.099127e-08 3.099127e-08 
# d01A_01B     d00B_00A     d11B_11A     d01B_01A 
# 3.099127e-08 3.099127e-08 3.099127e-08 3.099127e-08 

#####################################################################
# 9. Reconstruir rangos ancestrales (MarginReconGeoSSE)

# Estima el estado geográfico (y oculto) más probable en cada nodo interno del árbol.

recon_m_1 <- MarginReconGeoSSE(
        phy = phy.drop,
        data = regions.phy.df,
        f = f,
        pars = mejor_modelo_1$solution,
        hidden.states = hidden_geo,
        AIC = mejor_modelo_1$AIC
)
recon_m_3 <- MarginReconGeoSSE(
        phy = phy.drop,
        data = regions.phy.df,
        f = f,
        pars = mejor_modelo_3$solution,
        hidden.states = hidden_null,
        AIC = mejor_modelo_3$AIC
)

head(recon_m_3)
class(recon_m_3$node.mat)


recon_m_3$node.mat

#Identificar el estado más probable para cada nodo
j <- apply(recon_m_3$node.mat[, -1], 1, function(x){
        nombres <- colnames(recon_m_3$node.mat[, -1])
        nombres[which.max(x)]
})

prop.table(table(j))
# (00A)      (01A)      (11A) 
# 0.03703704 0.85185185 0.11111111 

# Rangos de transición
recon_m_3$phy$trans.rate
# (00) (11) (01)
# (00)   NA    0    1
# (11)    0   NA    2
# (01)    0    0   NA
recon_m_3$phy$Nnode
recon_m_3$phy$edge
recon_m_3$phy$edge.length

####################################################################
# 10 Visualizar tasas y estados en la filogenia

# Graficar la reconstrucción con tasas de diversificación neta
plot.geohisse.states(
        x = list(recon_m_3),
        rate.param = "net.div",
        # opciones: "turnover", "net.div", "speciation", "extinction","extinction.fraction"
        type = "fan",
        # "fan" o "phylogram"
        show.tip.label = TRUE,
        legend = TRUE,
        fsize = 0.8
)
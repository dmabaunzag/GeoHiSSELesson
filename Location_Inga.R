phy$tip.label

library(readxl)

# leer tabla de la distribución del género Inga
location.inga <- read_excel("Location_Inga.xlsx", sheet = "Table008 (Page 6-17)")
head(location.inga)
View(location.inga)

#resumen por región
regions <- table(location.inga$species, location.inga$region)
rownames(regions)

#Ajustar a los nombre del árbol filogenético
rownames(regions) <- sub("\r\n", "_", rownames(regions))

#Verificar si coinciden
matched.inga <- phy$tip.label %in% rownames(regions)

#Las que coinciden
phy$tip.label[which(matched.inga)]

# LAs que no
phy$tip.label[which(!matched.inga)]

sub("\r\n", "_", rownames(regions))[which(!sub("\r\n", "_", rownames(regions)) %in% phy$tip.label)]


regions.phy <- regions[which(matched.inga), c(1,3)]

rownames(regions.phy) %in% phy$tip.label



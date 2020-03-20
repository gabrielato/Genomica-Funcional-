## Alineamientos multiples en R 

#Librerias necesarias 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

###

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa") 

###

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

###

library(Biostrings)
library (msa)
library(seqinr)
library(ape)
library (ggtree)
##

#1. Secuencias concatenadas (en formato fasta)
#Si es una secuencia de DNA usamos readDNAStringSet 
#Si es Aminoacidos usamos readAAStringSet

readAAStringSet("K:/BIOINFO/insulins.fasta") -> insulinas


#2. Alineamiento de secuencias
#Hay tres formas: msa, msaMuscle y 
a_insulinas_muscle <- msaMuscle(insulinas) 
print (a_insulinas_muscle, show="alignment") #ver cómo quedo alineado

#3. Árbol 
#Se deben hacer tres pasos para tener un arbo: Convertirlo, tener los valores de las distancias y convertirlo en matriz 
a_convertido <- msaConvert(a_insulinas_muscle, type = "seqinr::alignment")
a_distancia <- dist.alignment (a_convertido, "identity") 
a_matriz <- as.matrix(a_distancia)

#Podemos hacer el árbol podemos hacerlo con el paquete "ape"
insulinatree_nj <- nj (a_matriz)
plot(insulinatree_nj)

#El arbol tambien se puede hacer con el paquete "ggtree"

ggtree(insulinatree_nj)

ggtree (insulinatree_nj, branch.length='none', color= "white") + 
  geom_nodepoint(color="orchid4", shape= 10, size=10) +
  geom_tippoint(color="royalblue3", shape=9, size=3) +
  theme_tree("tan")

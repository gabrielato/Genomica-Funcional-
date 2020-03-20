############# BIOINFORMATICA ###############
#### PRIMER PARCIAL #### 
getwd() #equivalente a pwd 
# esto unir la terminal con R y nos dice donde estamos
setwd("~/media/uaq/GAB/clase_bioinfo") #equivalente a cd
#cambiarme de carperta pero esto se hace en consola
list.files() #equivalente a ls -al
#Tambien se pueden usar comodines
list.files (pattern = ".R")
#Nos busca todos los que terminan en .R, se puede hacer para PDF, etc.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install() #hacerlo en windowns 


## Clase 19 de Ago 2019
#### Expresiones y asignaciones#### 
#la asignacion se queda de manera permanente por eso se pone <-
x <- rnorm(10)
x
#### Vectores ####
triste <- c() #vector vacio 
triste

feliz <- c(seq(1,10,1))
feliz 

feliz [-length(feliz)] #quita el ultimo elemento

hola <- c (1,2,3,4)
hola
bye <- rep (hola, each=2) #va a repetir los digitos del vector dos veces c/u
bye
hola [1:5] #pone los elementos de un vector desde el 1 al 5


#ejercicio de vectores
edades<-c(35,35,70,17,14)
nombres <-c("Jerry","Beth","Rick", "Summer","Morty")
names(edades)<-nombres
mean (edades [-2]) #el promedio de la edad sin Beth 

nuevoedad <- edades [-length(edades)] 
nuevoedad #quita a Morty del vector y creamos un nuevo vector a partir de esa opcion 
sort (nuevoedad) #ordena de menor a mayor

edades >= 75 #alguna edad es mayor a 75 años

edades <= 12 #alguna edad es menor a 12 años

edades >=12

edades <=20

any (edades >=12 & edades <=20) #buscar alguno que este entre 12 a 20

#Ejemplos biologicos 
genomeSize<-c(3234.83,2716.97,143.73,0.014281,12.1) #agrupa los distintos tamaños de genomas en un vector 
genomeSize*1e6 #para ver el tamaño con el tamaño mega (1 millon)

organismo<-c("Human","Mouse","Fruit Fly","Roundworm","Yeast")
# hacer el vector con el nombre de los organismos correspondientes al anterior genoma

names(genomeSize) <- organismo #estamos agrupando el vector genomas y el de organismos en un solo
genomeSize


#Ejercicio clase 2
edadesmicro <- c (22, 20, 23, 20, 19, 21, 20, 22, 24, 18)
nombresmicro <- c ("Gissel", "Carla", "Cesar", "Iliana", "Fer", "Mariana", "Dhamaris", "Jennifer", "Alejandra", "Brenda")
names (edadesmicro) <- nombresmicro
edadesmicro
max (edadesmicro)
min (edadesmicro)
mean (edadesmicro)
median (edadesmicro)
sd (edadesmicro)
length(edadesmicro)
edades_impares <-edadesmicro[c(1,3,5,7,9)]
edades_impares

#quitar el max y el min, asi que primero ordeno del menor al mayor
edadesmicro_1 <- sort(edadesmicro)
edadesmicro_1

hist (edadesmicro_1 [edadesmicro_1 !=min(edadesmicro_1)])
hist (edadesmicro_1 [edadesmicro_1 !=max(edadesmicro_1)])

#### Como ver si son iguales
x <- c(1,2,4)
y <- c(3,1,3)
any (x==y) #nos pone false porque no son iguales
any (x !=y) #nos pone true porque no son iguales
sum (x==y) #como no son iguales nos regresa un vector de 0

x1 <- seq (1,1000,1)
y2 <- seq (1,1000,1)
any (x1==y2)
sum (x1==y2)
#Va a regresar un vector  de 10000 1, porque todos son verdaderos

## Factores
factor_x<- c(1,2,2,3,1,2,3,3,1,2,3,3,1)
factor_x
as.factor(factor_x)
#Los niveles son los factores que exiten, entonces solo hay 1, 2,3
#ejercicio 

virus <- c("abutilon mosaic bolivia virus", "actinomyces virus", "alphapapilloma virus", "apple dimple fruit viroid", "autographa californica multiple nucleopolyhedrovirus", "bacillus phage phrodo", "banana streak my virus", "beihai hepe like virus", "ailuropoda melanoleuca papillomavirus", "aeropyrum pernix spindle shaped virus")

numero <-  c ( "NC_015045", "NC_009643", "NC_001357", "NC_003463", "NC_001623", "NC_031100", "NC_006955", "NC_032451", "NC_035201", "NC_028268")

setwd("~/GENOMAS_VIRUS")
list.files (pattern = ".fna")
virus <- list.files (pattern = ".fna")
#pone en un vector todos los genomas de los virus 


####### Biostrings ######


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

####

setwd("~/GENOMAS_VIRUS")
list.files()
library(Biostrings)
x <- readDNAStringSet("NC_015045.fna")
x

#Te dice informacion 
width(x) #tamaño
names (x) #nombre del organismo 
alphabetFrequency (x) # te dice la frencuencia que hay de cada numero de bases
translate (x)  #proteinas
#crear una secuencia 
DNAStringSet("ATGGT") -> y
y
names (y)
#Aqui se concateno varios organismos entonces
#todas estas opciones nos permiten analizar el 
#conjunto de archivos 
todos <- readDNAStringSet("todos.fna")
todos

names (todos)

############  MATRICES ##################
matriz1 <- matrix
(c(1,5,9,8,-4, 7),nrow=3,ncol=2) #como crear una matriz
# 1, 5. 9 va en la columna, uno abajo del otro 
matriz1

colnames(matriz1)<-LETTERS[1:2] #le cambia el nombre a las columnas
rownames(matriz1)<-letters[3:5] #le cambia el nombre a las filas

matriz1

matriz1 [,2] #solo muestra la columna 2

matriz1[,2]>3 #te dice si algunos componentes de la columna 2 es mayor a 3

### ejercicios 
# hacer una matriz con los siguientes datos 
micromatrix <- matrix (nrow = 5, ncol = 4)
micromatrix #primero creamos una matriz vacia para poder asignarle datos 

nombres <- c("Yanet", "Dulce", "Nataly", "Gaby", "Ale")
rownames (micromatrix) <- nombres 

columnas <- c("mes_de_cumple", "dia_de_cumple", "materias", "numero fav")
colnames(micromatrix) <- columnas

#Para cambiar la primera columna: mes de cumpleaños 
micromatrix [1,1] <- 8
micromatrix [2,1] <- 1
micromatrix [3,1] <- 3
micromatrix [4,1] <- 10
micromatrix [5,1] <- 12

#Para cambiar la segunda columna: dia de cumple 
micromatrix [1,2] <- 5
micromatrix [2,2] <- 12
micromatrix [3,2] <- 19
micromatrix [4,2] <- 15
micromatrix [5,2] <- 8

#para cambiar la tercera columna: materias
micromatrix [1,3] <- 5
micromatrix [2,3] <- 4
micromatrix [3,3] <- 4
micromatrix [4,3] <- 6
micromatrix [5,3] <- 5

#para cuarta columna: numero fav
micromatrix [1,4] <- 5
micromatrix [2,4] <- 3
micromatrix [3,4] <- 7
micromatrix [4,4] <- 15
micromatrix [5,4] <- 7

# calcular promedio de materias que llevan 
mean (micromatrix [,3])
#fechas que mas se repiten
table (micromatrix [,1])

#desviacion estandar de numero fav 
sd (micromatrix[,4])

##### SEGUNDO PARCIAL #####
#### Biostrings ####

#Para una unica secuencia 
dna1 <- DNAString("ACGT-N")
dna1 #Nos dice que es y nos pone que es una secuencia y cuantos nucleotidos tiene

#mas de 1 secuencia
dna2 <- DNAStringSet(c("ACGT", "GTCA", "GCTA"))
dna2

#- (para insercion) 
# N para cualquiera de los 4 nucleotidos

IUPAC_CODE_MAP # codificacion de todos los simbolos de IUPAC

dna2 [2] #nos selecciona la secuencia dos
dna2 [[2]][2:3] # de la secuencia 2 nos selecciona los nucleotidos 2 y 3

rev(dna2) #nos da las reversas de cada secuencia 

readDNAStringSet("E:/BIOINFO/secuencias_5.fasta") -> secuencias #lee el documento que hicimos en FASTA
secuencias # viene con los tamaños, los nombres, numero de acceso

names(secuencias)
alphabetFrequency(secuencias) #nos dice que tantos aminoacidos aparecen en la secuencia
# subseq  -> puede remplazar subsecuencias de inicio, fin de las secuencias o en un vector
subseq(secuencias)
sort (secuencias) #ordena alfabeticamente. 

#### Alineamiento de secuencias #### 
#cargar una taba
y <-read.table ("L:/ProteinTable167_161521.txt")
y #esto no se ve 

#descargar el formato tabular de la especie de nuestro interes
#inducada la longitud en aa 
Ecolix <-read.csv("E:/ProteinTable167_161521.txt", sep= "\t")
Ecolix
View(Ecolix)

length(Ecolix) #Para saber cuantas columnas tiene 

class (Ecolix) #para saber que es, en este caso es una data frame 

Ecolix$Stop #te da un vector que que incluye todos los stops de la tabla

#Ctrl 1 es para quedarnos aqui, en mi script
#Ctrl 2 es para pasarme a la consola 

sum(Ecolix$Strand=="+") #sumatoria de todos los verdaderos, estops representran los +, por lo cual existen 2063 + y los demas son los - 
sort (Ecolix$Length) #para ordenar del menos al mayor y ver cual es el más grande
#hay que seleccionar el ultimo que representa el más pequeño
sort (Ecolix$Length, decreasing = TRUE)
sort (Ecolix$Length, decreasing = TRUE)[1] ->grande
grande
max(Ecolix$Length) ->gengran #otra manera de saber cuanto mide el gen más grande
#ahora hay que saber quien es el más grande 
which (grande==Ecolix$Length) #el renglon en donde esta el más grande 
#cual de todas las longitudes representa a este gen mas grande
Ecolix [1937, "Locus"]
paste("longitud del gen mas grande", grande)
#como encontrar el gen mas grande y como se llama 
paste ("El gen más grande mide", grande, "aa" , "y se llama", locus)

##para el gen mas chico 
sort (Ecolix$Length) [1] ->chico 
wich (chico==Ecolix$Length)
Ecolix [1454, "Locus"] -> nombre
paste ("El gen más chicho se llama ", nombre, "y mide", chico, "aa")

## ejercicios extra 
median (Ecolix$Length) -> mediana #optener la mediana de los genes 
#ver si existe un gen que mida justo lo de la mediana 
sum (Ecolix$Length == "mediana")


#### lenguaje de seleccion ####
#valores booleanos, solo TRUE o FALSE 
#### CONDICIONAL #### 
x <- 5
if (x > 3) {print ("x es mayor que 3")}
#Lo imprime porque X es mayor a 3

#mas de 1 condicion
if (x < 3) {print ("x es mayor que 3")} else {
  print ("x es menor que 3")
}
# lo imprime porque X no es menor a 3 

#mas de 2 condiciones
if (x == 3) {print ("x es igual que 3")} else if (x < 3) {print ("x es menor que 3")} else if (x > 3) {print ("x es mayor que 3")}
#imprime la ultima porque ninguna es verdadera excepto la ultima 

#multiples condiciones
if (x >0) {
  x*5 
  print ("x es mayor que 0")
  print ("Hola")
}

if (x > 0) {
  print ("x es mayor que 0")
} else {print ("x es mayor que 0")}
#como x si es mayor que 0 solo imprime eso y no la segunda condicion 

ifelse(x>0, "x es mayor que 0", "x es menor que 0")
#en la primera ponemos la condicion, YES, NO 
#y como x es mayor que 0 imprime eso 

y <- -3
if (y>0) {
  print ("positivo")
} else if (y<0) {
  print ("negativo")
} else { 
  print("Cero")
}
#como es negativo unicamente imprime esa condicion 

dia <- 30
mes <- 1
if (dia>20 & mes <=2) {
  print ("Eres acuario")
}

# programa que te diga tu signo zodiacal 
dia <- 15
mes <- 10
if (dia > 24 & mes <=9) {
  print ("Eres libra")
} else if (dia < 22 & mes>=10) {
  print ("Eres libra")
}

# otra opcion 
if ((dia >= 24 & mes== 9 ) | (dia >= 22 & mes == 10)) {
  print ("eres libra")
} else { 
  print  ("no eres libra")
}
## con todos los signos zodiacales
dia <- 15
mes <- 8
if (dia > 21 & mes ==3) {
  print ("Eres aries")
} else if (dia < 20 & mes==4) {
  print ("Eres aries")
} else if (dia > 21 & mes ==4) {
  print ("Eres tauro")
} else if (dia < 20 & mes==5) {
  print ("Eres tauro")
} else if (dia > 21 & mes ==5) {
  print ("Eres Geminis")
} else if (dia < 21 & mes==6) {
  print ("Eres Geminis")
} else if (dia > 22 & mes ==6) {
  print ("Eres Cancer")
} else if (dia < 22 & mes ==7) {
  print ("Eres Cancer")
} else if (dia > 23 & mes ==7) {
  print ("Eres Leo")
} else if (dia < 23 & mes == 8) {
  print ("Eres Leo")
} else if (dia > 24 & mes ==8) {
  print ("Eres Virgo")
} else if (dia < 23 & mes ==9) {
  print ("Eres Virgo")
} else if (dia > 24 & mes == 9) {
  print ("Eres libra")
} else if (dia < 22 & mes == 10) {
  print ("Eres libra")
} else if (dia > 23 & mes == 10) {
  print ("Eres Escorpion")
} else if (dia < 22 & mes == 11) {
  print ("Eres Escorpion")
} else if (dia > 23 & mes == 11) {
  print ("Eres Sagitario")
} else if (dia < 21 & mes == 12) {
  print ("Eres Sagitario")
} else if (dia > 21 & mes == 12) {
  print ("Eres capricornio")
} else if (dia < 19 & mes == 1) {
  print ("Eres capricornio")
} else if (dia > 20 & mes == 1) {
  print ("Eres Acuario")
} else if (dia < 19 & mes == 2) {
  print ("Eres Acuario")
} else if (dia < 20 & mes == 2) {
  print ("Eres picis")
} else if (dia > 20 & mes == 2) {
  print ("Eres picis")
}

#### Multiple Sequence Alignment #### 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa") 

###
#Repasito de biostrings
library(Biostrings)
seq1 <- DNAString("AGTC") 
seq2 <- DNAString("ACTG") 
seq3 <- DNAString("ACCT")

readDNAStringSet("E:/BIOINFO/secuencias_5.fasta") -> seq4
alphabetFrequency(seq4) #Frecuencia de nucleotidos
dinucleotideFrequency(seq4) #Frecuencia de dinucleotidos 
reverseComplement(seq4)
reverse(seq4)
complement(seq4)

## alineamiento con secuencias pre cargadas
alineamiento1 <- system.file("examples","exampleAA.fasta",package="msa")
alineamiento1 <- readAAStringSet(alineamiento1)
primeralineamiento <- msa (alineamiento1)
print (primeralineamiento, show= "complete")
readDNAStringSet ("E:/BIOINFO/GENOMAS_VIRUS/NC_001357.fna", format="fasta") -> seqvirusAA
translate(seqvirusAA) -> seqAA
alphabetFrequency(seqAA)


## alineamiento diferente 

secuencia <-readAAStringSet(system.file("E:/BIOINFO/secuencias_5.fasta",package="msa"))
secuencia <-msa(secuencia) ## No corre porque tal vez vienen NA´s


#### Calcular probabilidad ####
genoma <- readDNAStringSet("K:/BIOINFO/GENOMAS_VIRUS/NC_001357.fna") 
frecuencia1 <- dinucleotideFrequency(genoma, as.prob=TRUE)
frecuencia1


matriz <- matrix(frecuencia1, nrow=4) #A las frecuencias que ya hicimos, sacar una matriz
matriz

frecuencia2 <- trinucleotideFrequency(genoma, as.prob = TRUE)
frecuencia2

##Para buscar una secuencia en el genoma, cuantas veces a aparece 
sequencia <- DNAString("ATGAG")
vmatchPattern(sequencia, genoma) # Te dice cuantas veces aparece esa secuencia en el genoma que aparece

##un no se que 
mean (tamanos) -> promedio 
which (tamanos > mean (tamanos)) -> indices

#### Estructuras de repeticion ####
#hay muchas, se manejan por ciclos 
#Ciclo, loops
#### CICLO FOR #### 
#Ciclo for: un numero de operaciones, un numero determinado de veces. Uno sabe cuantas veces debe hacer algo

for (i in 1:10)
{
  print ("Hola")
  print (i)
}

#Hacer un programa que diga Hola las veces de los nombres 
nombres <- c ("Isaac", "Brenda", "Perla", "Amy", "Ale", "Gaby")
nombres

for (i in 1:6)
{ 
  print (paste("Hola", nombres[i]))
  print (i)
}

#para cuando no queremos contar a mano 
for (i in 1:length(nombres))
{ 
  print (paste("Hola", nombres[i]))
  print (i)
}

#que sume todos los numeros del 1 al 10 
suma <- 0
for (i in 1:10)
{
  suma <- suma + i 
}
suma # vale 5050 la suma total 

colores <- c ("verde", "azul", "rojo", "morado", "azul", "dorado")
## ahora con colores fav 
for (i in 1:length(nombres))
{ 
  print (paste("Hola", nombres [i], "tu color fav es", colores[i]))
  print (i)
}

#"  ahora con secuencias 
readDNAStringSet("K:/BIOINFO/GENOMAS_VIRUS/NC_001357.fna") -> seq1
readDNAStringSet("K:/BIOINFO/GENOMAS_VIRUS/NC_001623.fna") -> seq2
readDNAStringSet("K:/BIOINFO/GENOMAS_VIRUS/NC_003463.fna") -> seq3
readDNAStringSet("K:/BIOINFO/GENOMAS_VIRUS/NC_006955.fna") -> seq4
readDNAStringSet("K:/BIOINFO/GENOMAS_VIRUS/NC_009643.fna") -> seq5

secuencias <- c (seq1, seq2, seq2, seq3, seq4, seq5)

names (secuencias) -> secuencias1
secuencias1

width(secuencias) -> secuencias2
secuencias2


for (i in 1: length(secuencias))
{
  print (paste ("El nombre de la secuencia", secuencias1 [i], "y mide nucleotidos", secuencias2 [i]) )
  print (i)
}

## ejercicio del examen 2: problema 2
sec_1 <- DNAString("AACTG")
sec_2 <- DNAString("AACTT")

diferencias <- 0 

for (i in 1:length(sec_1))
{
  if (sec_1[i]!=sec_2[i]) {
    diferencias <- diferencias + 1
  }
}
diferencias # solo hay una diferencia

## con el ejercicio del examen 

readAAStringSet("E:/BIOINFO/SEGUNDO PARCIAL ROBERTO/todosaa.fasta") -> secc
secc

### Para que nos de el nombre de todas las secuencias que existen en el archivo secc y tambien nos de el tamaño de aa
names (secc) -> nombressec
nombressec ## necesitamos hacer este objeto que incluya todos los nombres 

width(secc) -> tamañosec
tamañosec ## también un objeto donde vienen todos los tamaños 

length(secc)

for ( i in 1:length(secc) )
{
  print (paste ("El nombre de la secuencia", nombressec [i], "y mide nucleotidos", tamañosec [i]) )
  print (i)
  
}

### ahora hay que ver el numero de diferencias que existen entre todas las secuencias de secc que miden lo mismo o sea 110 aa (insulinas)
readDNAStringSet("E:/BIOINFO/insulinas/sequence(1).fasta") -> seq1
readDNAStringSet("E:/BIOINFO/insulinas/sequence(2).fasta") -> seq2

todos_seq <- c(seq1, seq2)
width(todos_seq)

diferencias <- 0 

for (i in 1:length(seq1))
{
  if (seq1[i]!=seq2[i]) {
    diferencias <- diferencias + 1
  }
}

diferencias ## esta secuencia difiere en 1 aa 


### utilizar el ciclo for para el ejercicio 5 de la tarea

poderes <- c (4,9,7,8.2,10,7.5,8,3,8.5, 3, 9)
superheroes <- c ("ironman", "capitanamarvel", "capitanamerica", "thor", "wanda", "thor", "wanda", "hulk", "drstrange", "blackwidow", "thanos")


for (i in 1:length(superheroes))
{ print (paste ("El superheroe", superheroes[i], "tiene un poder de", poderes [i]) )
  print (i)
  
}

### como lo hizo Roberto 
personajes <- c("Hulk", "Wanda", "Spiderman", "Hank")
poder <- c ("es verde","habilidades magicas","araña","nada")
valores <- c (7, 10, 8, 0)
super_heroes <- data.frame(personajes, poder, valores)
super_heroes

for (i in 1:length(poder)) {
  if (valores[i] >=8) {
    print (paste(personajes[i], "con el super poder", poder [i], "conlleva una gran bla bla"))
  }
}

#### GGTREE #### 
setwd("K:/BIOINFO/GENOMAS_VIRUS/")
secuencias <- list.files() #da todos los nombres de los atrchivos que estan en esa la carpeta 
secuencias
sec_reales <- readDNAStringSet(secuencias) # ahora lee todos los documentos dna
sec_reales

#LIBRERIAS
##
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
## 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("treeio")
## 
install.packages("ggimage" )
library(ggimage) 
install.packages("ggstance")
library(ggstance)
install.packages("ggnetwork")
library(ggnewscale)

##
## jugando con arboles
arbol1 <- system.file("extdata", "sample.nwk", package="treeio")
arbol <- read.tree(arbol1)
arbol
ggtree(arbol, ladderize=TRUE, branch.length="none") 
#el true es para voltear la imagen
ggtree (arbol, layout="slanted")

ggtree (arbol, branch.length='none', layout='circular', color= "white") + 
  geom_nodepoint(color="skyblue1", alpha=1/4, size=10) +
  geom_tippoint(color="wheat", shape=8, size=3) +
  theme_tree("black")

### AHORA CON SECUENCIAS REALES

readAAStringSet("K:/BIOINFO/insulins.fasta") -> insulinas
insulinas
library (msa)
library(seqinr)
library(ape)

### Aqui hacemos el arbol con ape 
a_insulinas_muscle <- msaMuscle(insulinas) #alineamiento

a_convertido <- msaConvert(a_insulinas_muscle, type = "seqinr::alignment")
a_distancia <- dist.alignment (a_convertido, "identity") 

a_matriz <- as.matrix(a_distancia)
insulinatree_nj <- nj (a_matriz)
plot(insulinatree_nj)

### gg tree
ggtree(insulinatree_nj)

ggtree (insulinatree_nj, branch.length='none', color= "white") + 
  geom_nodepoint(color="orchid4", shape= 10, size=10) +
  geom_tippoint(color="royalblue3", shape=9, size=3) +
  theme_tree("tan")

colors() ## para ver la lista de colores

#### CICLO WHILE ####
# se ejecuta algo solo si la condicion es cierta 
#hay que ponerle algo que determine la condicion 

x <- 3
while (x <= 5)
{ print (paste (x, " es menor que 5"))
  x <- x + 1
}
# como le pusimos un nuevio valor a x, va a repetir el ciclo hasta que el valor valga 5
#y ya no podrá cumplir la condicion
#se ejecuta tantas veces la condicion sea verdadera

####  interactuar con el usuario #### 
nombre <- readline(prompt = "Dime tu nombre: ")
signo <- readline (prompt = "dime tu signo: ")

print(paste ("Hola, ", nombre, "tu signo es", signo))

nombre <- "anonimo"
while (nombre !="NO") {
  nombre <- readline(prompt = "Dime tu nombre: ")
  signo <- readline (prompt = "dime tu signo: ")
  
  print(paste ("Hola, ", nombre, "tu signo es", signo))
}

### ejercicio que el limite de edad sea 100 años 
edad <- 0
edad <- as.numeric(edad)
while (edad !=100) {
  nombre <- readline (prompt = "Dime tu nombre: ")
  edad <- readline (prompt = "dime tu edad: ")
  
  print (paste ("Hola, ", nombre, "tu edad es", edad))
}

# solo va a parar cuando la persona de como edad 100. ahi se termina el ciclo 

####
#### EJERCICIOS PARA GANAR PUNTOS ####
####
## CICLOS FOR 
# tome 10 archivos de anotacion gff y dibuje un boxplot del tamaño de los genes para cada organismo 

setwd("E:/genoma/")
list.files() -> secuencias

sec <- readDNAStringSet(secuencias)

read.csv ("GCF_000005845.2_ASM584v2_genomic.gff.gz", sep= "\t") -> x
read.csv ("GCF_000006745.1_ASM674v1_genomic.gff.gz", sep= "\t") -> w
read.csv ("GCF_000006765.1_ASM676v1_genomic.gff.gz", sep= "\t") -> e
read.csv ("GCF_000008525.1_ASM852v1_genomic.gff.gz", sep= "\t") -> t
read.csv ("GCF_000009065.1_ASM906v1_genomic.gff.gz", sep= "\t") -> i
read.csv ("GCF_000009205.2_ASM920v2_genomic.gff.gz", sep= "\t") -> o
read.csv ("GCF_000013425.1_ASM1342v1_genomic.gff.gz", sep= "\t") -> p
read.csv ("GCF_000195955.2_ASM19595v2_genomic.gff.gz", sep= "\t") -> u
read.csv ("GCF_000195995.1_ASM19599v1_genomic.gff.gz", sep= "\t") -> ñ
read.csv ("GCF_001022195.1_ASM102219v1_genomic.gff.gz", sep= "\t") -> b

####
año <- 0
año <- as.numeric(año)
while (año !=2020) {
  año <- readline (prompt = "Dime tu año de nacimiento: ")
  
  if ((año== 1936) |(año == 1948)| (año== 1960) | (año== 1972) | (año == 1984) | (año == 1996) | (año == 2008 )){
    print ("RATA")
  } else if ((año== 1937) |(año == 1949)| (año== 1961) | (año== 1973) | (año == 1985) | (año == 1997) | (año == 2009 )) {
    print ("BUEY")
  } else if ((año== 1938) |(año == 1950)| (año== 1962) | (año== 1974) | (año == 1986) | (año == 1998) | (año == 2010 )) {
    print ("TIGRE")
  } else if ((año== 1939) |(año == 1951)| (año== 1963) | (año== 1975) | (año == 1987) | (año == 1999) | (año == 2011 )) {
    print ("CONEJO")
  } else if ((año== 1940) |(año == 1952)| (año== 1961) | (año== 1976) | (año == 1988) | (año == 2000) | (año == 2012 )) {
    print ("DRAGON")
  } else if ((año== 1941) |(año == 1949)| (año== 1962) | (año== 1977) | (año == 1989) | (año == 2001) | (año == 2013 )) {
    print ("SERPIENTE")
  } else if ((año== 1942) |(año == 1950)| (año== 1963) | (año== 1978) | (año == 1990) | (año == 2002) | (año == 2014 )) {
    print ("CABALLO")
  } else if ((año== 1943) |(año == 1951)| (año== 1964) | (año== 1979) | (año == 1991) | (año == 2003) | (año == 2015 )) {
    print ("CABRA")
  } else if ((año== 1944) |(año == 1952)| (año== 1965) | (año== 1980) | (año == 1992) | (año == 2004) | (año == 2017 )) {
    print ("MONO")
  } else if ((año== 1945) |(año == 1953)| (año== 1966) | (año== 1981) | (año == 1993) | (año == 2005) | (año == 2018 )) {
    print ("GALLO")
  } else if ((año== 1946) |(año == 1955)| (año== 1967) | (año== 1982) | (año == 1994) | (año == 2006) | (año == 2019 )) {
    print ("PERRO")
  } 
  
  print (paste ("Hola ", año))
}

readDNAStringSet("E:/BIOINFO/GENOMAS_VIRUS/NC_001357.fna") -> dna
dna

rna <- gsub ("T", "U", dna)
rna

matchPattern("UAG", rna)
matchPattern("UGA", rna)
matchPattern("UAA", rna)
matchPattern("AUG", rna)

cuadrado <- 0

for (i in 1:100)
{
  cuadrado <-  (sqrt(cuadrado)) + i
}

cuadrado

####
####
####

##### ggtree ####
## paquetes a descargar 
library (ggtree)

install.packages("phytools")
library (phytools)

#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EBImage")
library (EBImage)
#
install.packages("ggimage")
library (ggimage)

# Exportar un texto que hicimos en mega, formato newick 
arbol_newick <- read.newick("K:/BIOINFO/insulinas_arbol_nwk")
# usamos la funcion de read.newich, el cual nos permite leer ese archivo de texto del árbol que alineamos en MEGA e hicimos un arbol
# este arbol de insulinas, lo hicimos con alineamiento de clustal 
# y el arbol con maxima parsimonia 
plot(arbol_newick)

arbol_rec <- ggtree(arbol_newick, layout="rectangular")
plot (arbol_rec)

arbol_sla <- ggtree (arbol_newick, layout="slanted")
plot (arbol_sla)

arbol_cir <- ggtree (arbol_newick, layout = "circular")
plot (arbol_cir)

arbol_rad <- ggtree (arbol_newick, layout = "radial")
plot (arbol_rad)

arbol_cir <- ggtree (arbol_newick, layout = "circular")
plot (arbol_cir)

arbolsraiz <- ggtree (arbol_newick, layout = "unrooted")
plot (arbolsraiz)

multiplot(arbol_cir, arbol_rad, arbol_rec, arbol_sla, arbolsraiz, ncol= 2) 
# para visualizar todos los arboles juntos 


arbol_insulinas <- ggtree(arbol_newick) +
  geom_label2(aes(subset=!isTip, label=node)) + 
  geom_tiplab (size=2,  color = "black")
plot (arbol_insulinas) 


arbol_colores <- arbol_insulinas+ 
  xlim (0, 1) +  
  geom_hilight(node= 19, fill= "turquoise1", alpha = 0.5) + 
  geom_hilight(node= 22, fill= "green", alpha = 0.5) + 
  geom_hilight(node= 16, fill = "pink", alpha = 0.3) +
  geom_hilight(node= 25, fill = "tomato", alpha = 0.5) +
  geom_hilight(node= 24, fill = "turquoise1", alpha = 0.5) 

plot (arbol_colores)

##### Funciones definidas por el usuario #####
# asignamos funciones a objetos 

#nombre_funcion <- fuction (argumentos)
# {
# definir la funcion 
#}

##
mi_primera_funcion <- function(n)
{ for (i in 1:n)
{print (i)
}
}

mi_primera_funcion(10)
# va a correr los 10 numeors del 1 al 10 como esta en la funcion 

##
mi_segunda_funcion <- function(inicio, fin)
{
  for (i in inicio:fin) {
    print(i)
  }
}  

mi_segunda_funcion(10,30)
# va a correr de 10 al 30, ya que fue definido dos veces

## 
#Ejercicio: calcular la distancia de Halmin: tomar el numero de diferencias de dos secuencias, divididas por el numero de nucleotidos
sec_1 <- DNAString("AACTG")
sec_2 <- DNAString("AACTT")

diferencias <- 0

distancia <- function (sec_1, sec_2)
{diferencias <- 0
for (i in 1:length(sec_1)) {
  if (sec_1[i]!=sec_2[i]) {
    diferencias <- diferencias + 1}} 
return (diferencias)

}

distancia (sec_1, sec_2) # solo hay una diferencia

sec_3 <- DNAString ("AGTAT")
sec_4 <- DNAString("AACTT")

#otra funcion que divide el numero de diferencias entre el numero de nucleotidos totales
distancia_2 <- function (sec_1, sec_2)
{diferencias <- 0
{for (i in 1:length(sec_1)) {
  if (sec_1[i]!=sec_2[i]) {
    diferencias <- diferencias + 1}}} 
print (paste (diferencias, "de", length(sec_1)))
}

distancia_2(sec_3, sec_4) 

# 1) a partir de dos secuencias que generen un score de alineamiento 
# 2) dos nombres de mis compañeras y calcular la diferencia

# 1) 
sec_3 <- DNAString ("AGTAT")
sec_4 <- DNAString("AACTT")

alineamiento <- function(sec_1, sec_2) {
  return (pairwiseAlignment(sec_1, sec_2)) 
} 

alineamiento(sec_3, sec_4)

# 2) 
#para nombres con el mismo tamaño 

diferencias <- function (x, y) {
  diferencias <- 0 
  for (i in 1:length (y)){
    if (x[i] != y[i]) {
      diferencias <- diferencias + 1 
    }
  }
  return (diferencias)
}

nom1 <- AAString("Carla")
nom2 <- AAString ("Perla")

diferencias (nom1, nom2)

# para nombres de diferente tamaño 
nom3 <- AAString ("Alejandra")
nom4 <- AAString("Gabriela")



diferencias1 <- function (x, y) {
  diferencias <- 0 
  minimo <- min (length(x), length (y))
  for (i in 1:minimo){
    if (x[i] != y[i]) {
      diferencias <- diferencias + 1 
    }
  }
  return (diferencias)
}

diferencias1 (nom3, nom4)



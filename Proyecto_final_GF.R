####### Proyecto final genomica ########
########################################
#Gabriela Torres & Dhamaris Quevedo ####
########################################

## LIBRERIAS QUE SE VAN A NECESITAR
library(dada2)

## PASO 1 #### 
#Necesitamos subir las secuencias FASTQ que seleccionamos y descargamos previamente

setwd ("E:/PROYECTO FINAL/") #cambiamos al directorio en donde se encuentra la carpeta de las secuencias
path <-"SEQ/" #se asigna a un objeto, la carpeta en donde estan las secuencias 

list.files(path) #comprobamos que las secuencias se encuentren en el objeto path 

#Como se descargaron dos secuencias por muestra, es decir, la secuencia que corresponden a los
#al Forward (terminacion _1) y los reverse (terminacion _2), con la intencion de separar en dos
#objetos distintos, las secuencias correspondientes a cada tipo de muestra.
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

## PASO 2 ####
#Observar la calidad de las lecturas
#Esto se hace para saber si necesitamos recortar o filtrar las secuencias
#se hace tanto para las secuencias F y R

plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])

## PASO 3 ####
#Recortar secuencias 

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

#Primero hay que asignarles el nombre a cada muestra, el nombre filtrado
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

filtFs
filtRs #comprobamos que tengan el nombre adecuado 
#y lo siguiente a realizar es establecer los parametros para recortar la secuencia
# Aqui se establece el n�mero m�ximo de errores permitidos en una lectura

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
#La funcion filterandtrim se encarga de remover los  primers 
#y recortar / filtrar las sec para la calidad
head(out) #comprobar que se haya hecho bien

## PASO 4 ####
#Tasas de error

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

#Esta funcion se encarga de generar un modelo de error de nuestros datos

save(errF, file = "rF.RData")
save(errR, file = "rR.RData")
#Nota: como son funciones que tomaron muchisimo tiempo en correr, se procedio a 
#guardarlas como una base de datos de R, para poder subirla despu�s sin la necesidad de 
#volver a correrlo. 

#Hay que visualizae estas tasas de error
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


## PASO 5 ####
#Inferencia de muestra 

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#esta funci�n es para inferir las variantes de secuencia de amplic�n  en lecturas
#directas (F) e inversas (R) de forma independiente


save(dadaFs, file = "dFs.RData")
save(dadaRs, file = "dRs.RData")
#Guardamos el objeto, ya que tardo en cargar

#despues hay que ver que infirio el aloritmo para ambas secuencias
dadaFs[[1]]
dadaRs[[1]]

## PASO 6 ####
#Combinar las lecturas
#Esto se hace con la intencion de reducir el ruido
#La fusion solo se completa si se superponen en al menos 12 bases 
#y son id�nticas entre s� en la regi�n de superposici�n

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#Esta funcion se encarga de fusionar lecturas  
#para refinar a�n m�s las variantes de secuencia de amplic�n

head(mergers[[1]]) #este es un objeto que ess una lista de la base de datos de cada muestra

## PASO 7 ####
#Tabla de secuencia
#Esta es una tabla de las variantes de secuencia de amplic�n 
#y proporciona mayor resolucion que las OTUS

seqtab <- makeSequenceTable(mergers)
#Funcion que genera una tabla de conteo
dim(seqtab)

table(nchar(getSequences(seqtab))) #observar las variantes
#esta tabla es una matrix con filas (las muestras), y columnas (las variantes de secuencias). 


##PASO 8####
#Eliminar las quimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Esta funcion sirve para detectar y eliminar quimeras

dim(seqtab.nochim)
#vemos las dimensiones, para saber que porcentaje de quimeras es de acuerdo a las variantes

sum(seqtab.nochim)/sum(seqtab)
#Pero se debe hacer este paso, para realmente ver la abundancia de estas variantes

## PASO 9 ####
#Ver la cantidad de lecturas que quedaron despues de realizar todos los procesos

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

## PASO 10 ####
#Asignar taxonom�a 

#Comparaci�n de nuestras secuencias con un conjunto de secuencias de referencia de taxonom�a de la base de datos de Silva 

taxa <- assignTaxonomy(seqtab.nochim,refFasta = "~/genomica funcional/proyecto/silva_nr_v132_train_set.fa.gz",multithread=TRUE)

#con base en la coincidencia exacta de nuestras secuencias contra la base de datos de Silva, se genera 
#una clasificaci�n para determinar g�nero y especie 

taxa <- addSpecies(taxa, "~/genomica funcional/proyecto/silva_species_assignment_v132.fa.gz")


taxa.print <- taxa #Se eliminan los nombres de las secuencias para que sea m�s f�cil su visualizaci�n 
rownames(taxa.print) <- NULL
head(taxa.print) #visualiza nuestras secuencias, reino, phylum, clase, orden, familia, g�nero y especie 


save(taxa,file = "taxa.RData") ##### Guardamos este objeto como un documento de R ya que tard� demasiado en correr


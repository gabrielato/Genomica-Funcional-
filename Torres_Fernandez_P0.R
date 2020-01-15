#### PAQUETES ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

library(BiocVersion)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")

#### Tarea 1 #### 

# Problema 1
#traducir una secuencia RNA a AA
readRNAStringSet ("K:/GENOMICA FUNCIONAL/first.fasta", format="fasta") -> RNA #utilizando la 

translate(RNA) -> AA
#paquete especializado para traducir de RNA a Aminoacidos

#Problema 2
#Counting DNA Nucleotides
DNAString("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC")-> DNA
alphabetFrequency(DNA) # libereria especializada para el conteo

#Problema 3
#Computing GC Content
DNAString ("CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT")-> DNA1
dinucleotideFrequency(DNA1) -> frecuencia
frecuencia # vemos que se enuentra 5 veces 
length(DNA1) #vemos que el tamaño es de 87

(5/87) * 100 #porcentaje 

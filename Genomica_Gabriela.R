### GENOMICA FUNCIONAL ###

########### PRIMER PARCIAL ###################
#### PAQUETES ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

library(BiocVersion)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")


#### IGRAPH ####
install.packages("igraph")
library (igraph)
#   REDES  #### 

red <- make_empty_graph(n=5,directed = TRUE) # primero hay que hacer una vacia 
#V(red) es igual a todos losVertex(g)

V(red)$color ="red"  #para poner color
V(red)$shape ="sphere" #para figura

plot(red)

#ahora hay que darle conexiones,para eso se hace con un vector
red_1 <- add.edges(red, c(1,2, 1,3, 2,4, 3,4, 4,5, 1,4)) #esto son las conexiones
# las conexiones 1 conecta con 2 y asi 

plot (red_1)

#si despues quisieramos agregar otro vertice se hace

red_1 <- add.vertices(red_1, 1,  color= "blue", shape= "square") 
plot (red_1)

#y ahora hay que conectarlo con nuestra red 
red_1 <- add.edges(red_1, c(4, 6))
plot (red_1)

# para hacer un cambio de numeros a letras
V(red_1)$name <- LETTERS [1:6]
plot (red_1)

#### Degree ####
degree (red_1) #te pone con cuantos comparte de cada objeto

hist (degree(red_1),col = "pink") #graficar en histograma el degree

plot(red_1,layout=layout_nicely,vertex.size=degree(red_1,V(red_1),"in")*15+15,
     vertex.label.dist=0.5,edge.arrow.size=0.5) # lo grafica y entre más degree tenga, más grande la esfera

#head map del degree
matriz <-as.matrix(get.adjacency(red_1))
heatmap(matriz, Rowv=NA, Colv="Rowv")

#############################
## Tablas 
read.csv("k:/GENOMICA FUNCIONAL/genemarks.20200122.110338.1047.gms.out")

read.table("K:/GENOMICA FUNCIONAL/genemarks.20200122.110338.1047.gms.out", sep = "\t") -> tabla_virus

corona_virus <- readDNAStringSet("k:/GENOMICA FUNCIONAL/Wuhan_seafood_virus.fasta")
corona_virus

translate(corona_virus)
###############################


### RED: bajar una matriz y hacer tu red
read.csv("K:/GENOMICA FUNCIONAL/Red_de_Amigas_GF_2020 - Sheet1(1).csv") ->red_micro
red_micro

#hay que convertir los renglones con nombres

row.names(red_micro) <- red_micro[,1]

red_micro <- red_micro [,-1]
View (red_micro)

red_micro <- as.matrix(red_micro) #hay que transformarlo a matriz

red_friends <- graph_from_adjacency_matrix(red_micro, mode = "directed")
plot (red_friends, edge.arrow.size=0.3)

degree (red_friends, mode = "in") -> degree_amigos #que te consideran sus amigos
degree  (red_friends,mode= "out") -> degree_amigos_1 #que consideras sus amigos

mean (degree_amigos) -> promedio

# para ver grupos 
sp <- cluster_spinglass(red_friends)
plot(sp,red_friends, edge.arrow.size= 0.3)

#### Trayectoria y distancia ####
g<-random.graph.game(10,0.5) #para una red aleatoria de 10 nodos

plot(g)

distancia<-shortest.paths(g) #nos va a poner las distancias (los pasos que hay de un nodo a otro)
distancia #lo pone como una matriz 

heatmap(distancia)

distances(g) #lo mismo que poner shortest.paths
diameter(g) #diametro,la distancia MAS garnde, es 3

mean_distance(g) #el primedio de las distancias de la red g

#### Clustering Coefficient ####
transitivity(g) #es el promedio de la transitividad de TODOS los nodos. 

transitivity(g, type= "local") # Pone la transitividad de cada nodo. LOCAL

# ejercicios con la red de amigos 
plot (red_friends)

d1 <- shortest.paths(red_friends, mode= "in") #Consideras TUS amigos
d2 <- distances(red_friends, mode= "out") #Consideran SUS amigos 
heatmap(d1)
heatmap (d2)

mean_distance(red_friends)

shortest.paths(red_friends, "Gabriela", "Iliana") #la distancia que existe entre esos dos nodos
shortest.paths(red_friends, "Gabriela", "Omar") #como no estamos conectados, hay un nodo entre esos dos nodos

diameter(red_friends) #la mayor distancia es 3

transitivity(red_friends)

transitivity(red_friends, type= "local")


#### Free-Scale Property ####
barabasi.game() #comando para generar una red sin escala

sf <- barabasi.game(1000, directed = FALSE)
sf #ver las interacciones 
plot (sf,vertex.label=NA,edge.arrow.size=0.3) # para que no tenga los nombres

degree (sf) -> d_sf
hist (d_sf) #lamayoria de los nodos tienen 1 conexion 

boxplot (d_sf)

## ley de potencias ##
fit_sf <- fit_power_law(d_sf+1, 10) #distribucion de los datos es una ley de potencias
fit_sf$alpha #la potencia que calcula el ajuste. potencia ajustada a la distribucion 


#### ejercicios ####
g10<-barabasi.game(10,directed = FALSE)
g100<-barabasi.game(100,directed = FALSE)
g1K<-barabasi.game(1000,directed = FALSE)
g2<-random.graph.game(1000,0.20)
g3<-sample_smallworld(1,1000,p=0.2,nei=3)
g4<-make_graph("Zachary")

##
#funciones 
degree_hist <- function (n) {
  hist (degree (n))
}

degree_boxplot <- function (n) {
  boxplot (degree (n))
}

m_m_degree <- function(n) {
  paste ( "el promedio", mean (degree(n)), "y la media", median (degree(n)))
}

plot_dd <- function(n){
  plot (degree.distribution(g10))
}

##
g10
degree(g10) 
degree_hist(g10)
plot(degree.distribution(g10))
degree_boxplot(g10)
m_m_degree(g10)

g100
degree_hist(g100)
plot_dd(g100)
degree_boxplot(g100)
m_m_degree(g100)

g2
degree_hist(g2)
plot_dd(g2)
degree_boxplot(g2)
m_m_degree(g2)

g3
degree_hist(g3)
plot_dd(g3)
degree_boxplot(g3)
m_m_degree(g3)

g4
degree_hist(g4)
plot_dd(g4)
degree_boxplot(g4)
m_m_degree(g4)

######

##### Robustes de redes ####
g <- make_ring(10) %>%
  set_vertex_attr("name", value = LETTERS[1:10])

V(g) #paraver los nombres de los vertices de la red

g <- delete_vertices(g, c(1,5)) %>%
  delete_vertices("B") ## eliminar el vertice B 
g

plot (g)

#ejercicio y hay que quitar un numero al azar
g100<-barabasi.game(100,directed = FALSE)
g1K<-barabasi.game(1000,directed = FALSE)
g2<-random.graph.game(1000,0.20)
g3<-sample_smallworld(1,1000,p=0.2,nei=3)

sample(1:100,1) #que nos den un numero al azar del 1 al 100

g100<-delete.vertices(g100, sample(1:100,1))
plot(g100)
g100<-delete.vertices(g100, sample(1:100,1))
plot(g100)
g100<-delete.vertices(g100, sample(1:100,1))
plot(g100)
#y asi nos vamos hasta quitarle los 10 vertices, aunque sería mas facil con un ciclo for 


#### REDES ALEATORIAS ####
ra <- random.graph.game(1000,0.10)
degree(ra)
ra_d <-shortest.paths(ra)
mean(ra_d)
diameter(ra)

#ahora que nos quite
for (i in 1:100) 
{
  rad <- delete.vertices(g1000, sample(1:(1000-i),1))
  print (mean_distance(rad))
  print(diameter(rad))
}

############# SEGUNDO PARCIAL ######################


#### CENTRALIDAD ####
#para la red random
degree (ra) -> ra_degree
eccentricity(ra) -> ra_ec
closeness(ra) -> ra_clo
betweenness(ra) -> ra_bet

sort(ra_degree, decreasing = FALSE)
sort(ra_ec, decreasing = FALSE)
sort(ra_clo, decreasing = FALSE)
sort(ra_bet, decreasing = FALSE)

#para saber cuales
which(sort(ra_degree, decreasing = TRUE) [1]==ra_degree) #este es el numero del nodo mas conectado
sort(ra_degree, decreasing = TRUE) [1:10] #estos son los degrees de los 10 nodos mas conectados

## ahora con la red de amigos
read.csv("K:/GENOMICA FUNCIONAL/archivos/Red_de_Amigas_GF_2020 - Sheet1(1).csv") ->red_micro
row.names(red_micro) <- red_micro[,1]
red_micro <- red_micro [,-1]
red_micro <- as.matrix(red_micro) #hay que transformarlo a matriz
red_friends <- graph_from_adjacency_matrix(red_micro, mode = "directed")

degree(red_friends,mode= "in") ->bff_degree
degree(red_friends,mode= "out") ->bff_degree_out
eccentricity(red_friends, mode = "in") -> bff_ecc
eccentricity(red_friends, mode = "out") -> bff_ecc_out
closeness(red_friends, mode= "in") ->bff_clo
closeness(red_friends, mode= "out") ->bff_clo_out
betweenness(red_friends, mode="in") -> bff_bet
betweenness(red_friends, mode="out") -> bff_bet_out

sort (bff_degree, decreasing = TRUE)[1:10]
sort (bff_degree_out, decreasing = TRUE)[1:10]
sort (bff_ecc, decreasing = FALSE)[1:10]
sort (bff_clo, decreasing = TRUE)[1:10]

#### metodos de cluster ####
read.csv("~/GENOMICA_FUNCIONAL_GABRIELA/Red_de_Amigas_GF_2020 - Sheet1.csv") ->red_micro
row.names(red_micro) <- red_micro[,1]
red_micro <- red_micro [,-1]
red_micro <- as.matrix(red_micro) #hay que transformarlo a matriz
red_friends <- graph_from_adjacency_matrix(red_micro, mode = "directed")

#PARA REDES DURIGIDAS
ed <- cluster_edge_betweenness(red_friends, directed = TRUE)
ed
plot(ed, red_friends)

pg <- cluster_spinglass(red_friends)
pg
plot (pg, red_friends)

wt <- cluster_walktrap(red_friends)
wt
plot(wt,red_friends)

#### REDES PESADAS ####
#redes que tiene valores  
read.csv("~/GENOMICA_FUNCIONAL_GABRIELA/Weighted_Network_GF_2020 (Responses) - Form Responses 1.csv") -> redseries

row.names(redseries) <- redseries[,1]
redseries<- redseries[,-1]

View(redseries)

red_series <- as.matrix(redseries)
red_series

## perfil de gustos
plot (red_series["GABRIELA",], type="l")

cor (red_series["GABRIELA",], red_series["ILIANA",])
cor (red_series["ISAAC",], red_series["BRENDA",])

length(red_series[,1])

row.names(red_series) -> x
x
length(x)

for (i in 1:length(x))
{
  cor (red_series["GABRIELA",], red_series[i,]) ->r
  print(paste( "los valores de",x[i], "son",r))
} 

#matriz transpuesta para hacer un heat´map
cor (t(red_series)) -> perfiles
heatmap(perfiles)
View (perfiles)
#ahora hay que hacer esa matriz una red
perfiles <- perfiles -min(perfiles)
diag(perfiles) <-0

red <- graph_from_adjacency_matrix(perfiles,mode="undirected", weighted = TRUE)
plot(red)

#ahora en cluster
cluster_edge_betweenness(red, weights = E(red)$weight) ->red_ed
plot(red_ed, red)


#### REDES BOOLEANAS (COEXPRESION) ####
#Descargar en WGCNA y BoolNet
library(BoolNet)

loadNetwork("reg.txt") -> redV
redV #No corre, hacerlo a mano 

sink("regla.bn")
cat("targets, factors\n")
cat("V1, V1 | ! V3\n")
cat("V2, V1 & V3\n")
cat("V3, V2\n")
sink()

loadNetwork("regla.bn") ->reglas

atractores <- getAttractors(reglas)
atractores

plotAttractors(atractores)

redc <- data("cellcycle")

cellcycle

getAttractors(cellcycle) -> gt

plotAttractors(gt)


sink("operonlac.bn")
cat("targets, factors\n")
cat("M, !R\n")
cat("P, M\n")
cat("B, M\n")
cat("R, !A\n")
cat("A, B & L\n")
cat("L, P\n")
sink()

loadNetwork("operonlac.bn") -> operonlac

operonlac
plotNetworkWiring(operonlac) #no pone diferencias entre inhibir y no

getAttractors(operonlac) -> atractor
atractor
plotAttractors(atractor)
#esto demuestra 000100 (63 estados) y 111011 (1 estado)


## ahora queremos ver los estados en 000000 y 111111
stateTransition(operonlac, c(1,1,1,1,1,1)) #cuando todos estan prendidos se van a ese estado
stateTransition(operonlac, c(0,0,0,0,0,0)) #cuando estan apagados se van a este estado

plotStateGraph(atractor)


### REDES PESADAS DE CORRELACIONES#
library(WGCNA) #para redes muy grandes, hasta de 200000 genes

install.packages("BiocManager")
BiocManager::install("WGCNA")

#tutorial para estas redes
# cargar la base de datos
read.csv("K:/GENOMICA_FUNCIONAL_GABRIELA/bases de datos/ClinicalTraits.csv") -> ClinicalTraits
read.csv("K:/GENOMICA_FUNCIONAL_GABRIELA/bases de datos/GeneAnnotation.csv") -> Gene_Annotation
read.csv ("K:/GENOMICA_FUNCIONAL_GABRIELA/bases de datos/LiverFemale3600.csv") -> femData

library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Take a quick look at what is in the data set:
dim(femData)
names(femData)

datExpr0 = as.data.frame(t(femData[, -c(1:8)]))
names(datExpr0) = femData$substanceBXH
rownames(datExpr0) = names(femData)[-c(1:8)]
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

######## CYSTOCAPE #############
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RCy3")

library(RCy3)
cystocapePing() #para reconocer que esta conectado con cystocape


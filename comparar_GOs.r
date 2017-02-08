past_dir = getwd()
dir = "/media/julen/TOSHIBA_EXT/datos/paper/"
setwd(dir)

# Levadura
library("org.Sc.sgd.db")
library("gProfileR")
library("org.Hs.eg.db")
library("impute")
library("data.table")

############################ FUNCIONES #################################
# Preprocess like in SOTA scripts
preprocessData <- function(matriz) {
  matriz = matriz[,c("TR", "RS", "RA", "TLRi", "PS", "PA")]
  
  # Guardamos solo los genes de los que tengamos info en al menos 4 variables, 
  # y de distintas moleculas
  # Para SOTA utilizamos TLRi
  matriz1 = matriz[,c("TR","RS","RA")]
  posiciones= apply(!is.na(matriz1), 1, sum)
  pos = posiciones >= 2
  print("Genes con 2 de 3 variables en RNA")
  table(pos)
  matriz1 = matriz1[pos,]
  
  temporal = impute.knn(as.matrix(matriz1) ,k = 10, rowmax = 0.5, colmax = 0.8, 
                        maxp = 1500, rng.seed=362436069)
  matriz1 = as.data.frame(temporal$data)
  
  matriz2 = matriz[,c("TLRi","PS","PA")]
  posiciones= apply(!is.na(matriz2), 1, sum)
  pos = posiciones >= 2
  print("Genes con 2 de 3 variables en RNA")
  table(pos)
  matriz2 = matriz2[pos,]
  # Generamos valores faltantes con el knn
  temporal = impute.knn(as.matrix(matriz2) ,k = 10, rowmax = 0.5, colmax = 0.8, 
                        maxp = 1500, rng.seed=362436069)
  matriz2 = as.data.frame(temporal$data)
  
  # Juntamos las dos matrices
  pos = rownames(matriz1) %in% rownames(matriz2)
  table(pos)
  Matriz_6 = data.frame(row.names = rownames(matriz1)[pos])
  Matriz_6[,c("TR","RS","RA")] = matriz1[rownames(Matriz_6),c("TR","RS","RA")]
  Matriz_6[,c("TLRi","PS","PA")] = matriz2[rownames(Matriz_6),c("TLRi","PS","PA")]
  
  # Pasamos a matriz asegurandonos de que si alguna columna esta como factor, 
  #no se guarde su indice
  indx <- sapply(Matriz_6, is.factor)
  if (sum(indx) > 0) {
    Matriz_6[indx] <- lapply(Matriz_6[indx], 
                             function(x) as.numeric(as.character(x)))
  }
  Matriz_6 = as.matrix(Matriz_6)
  
  # Generamos ranking
  for(col in 1:ncol(Matriz_6)) {
    Matriz_6[, col] = rank(Matriz_6[, col], na.last = "keep")
  }
  # Lo pasamos a escala de 0 a 1
  Matriz_6 = t(t(Matriz_6)/apply(Matriz_6,2,max, na.rm = TRUE))
  return(Matriz_6)
}


# Preprocess keeping genes with data in two variables
preprocessData2 <- function(matriz) {
  # Guardamos solo los genes de los que tengamos info en al menos 2 variables
  Matriz_GO = matriz[,c("TR", "RS", "RA", "TLRi", "PS", "PA")]
  posiciones= apply(!is.na(Matriz_GO), 1, sum)
  pos = posiciones >= 2
  table(pos)
  Matriz_GO = Matriz_GO[pos,]
  
  # Generamos ranking
  Matriz_rank = Matriz_GO
  for(col in 1:ncol(Matriz_GO)) {
    Matriz_rank[, col] = rank(Matriz_GO[, col], na.last = "keep")
  }
  # Lo pasamos a escala de 0 a 1
  Matriz_rank = t(t(Matriz_rank)/apply(Matriz_rank,2,max, na.rm = TRUE))
  
  return(Matriz_rank)
}


# Function to obtain GO term from GO ID in HeLA and yeast
GoTermFinder <- function(common, GOtabl1, dataset, specie) {
  mart <- useMart(biomart = "ensembl", dataset = dataset)
  # Get associated term from BioMart
  goInfo <- getBM(attributes = c("go_id", "name_1006"), filters = "go_id",
                  values = common, mart = mart)
  # store in table
  GOtabl1[goInfo[,1],specie] <- goInfo[,2]
  return(GOtabl1)
}


# Function to obtain genes associated to a GO ID
GoGeneFinder <- function(common, dataset, id) {
  mart <- useMart(biomart = "ensembl", dataset = dataset)
  # Get associated term from BioMart
  goInfo <- getBM(attributes = c("go_id", id), filters = "go_id",
                  values = common, mart = mart)
  # store in table
  return(goInfo)
}

# Calculate standard error
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
############################# CODIGO #####################################
filtPercent <- 0.7
mini <- 50
maxi <- 250

################# Yeast ####################
gos <- org.Sc.sgdGO2ORF
gos <- as.list(gos)
# get genes on each GO
ngenes <- sapply(gos, function (x) length(x))
# Select only the ones with 50 or more and 250 or less
pos <- ngenes >= mini & ngenes <= maxi
gos <- gos[pos]
# Filter by presence of these genes in our dataset (say 70%)
lev <- readRDS(file="levadura/datos_cargar/matriz_lev.rds")
# lev <- preprocessData(lev)
lev <- preprocessData2(lev)
levGenes <- rownames(lev)
safe <- sapply(gos, function (x) {
          longi <- length(x)
          presence <- x %in% levGenes
          percent <- sum(presence) / longi
          if (percent >= filtPercent) {
            TRUE
          } else {
            FALSE
          }
        })

# Safe list of GOs
lgosYeast <- names(gos[safe])



############## HeLa ####################
# filtPercent <- 0.5

# Cambiamos los IDs a entrez
# org.Hs.egGO2ALLEGS is an R object that provides mappings between a given GO identifier and
# all of the Entrez Gene identifiers annotated at that GO term OR TO ONE OF IT’S CHILD NODES
# in the GO ontology. Thus, this mapping is much larger and more inclusive than org.Hs.egGO2EG.
gos <- org.Hs.egGO2EG
gos <- as.list(gos)
# get genes on each GO
ngenes <- sapply(gos, function (x) length(x))
# Select only the ones with 50 or more and 250 or less
pos <- ngenes >= mini & ngenes <= maxi
gos <- gos[pos]
# Open HeLa dataset
hela <- readRDS(file="HeLa/datos_cargar/matriz_HeLa.rds")
# hela <- preprocessData(hela)
hela <- preprocessData2(hela)
helaGenes <- rownames(hela)
# Convert symbol genes to entrez
xx1 <- as.list(org.Hs.egALIAS2EG)
listaGenes1 = xx1[helaGenes]
names(listaGenes1) = NULL
catego1 = unlist(listaGenes1)
helaGenes = unique(catego1)
# Filter by presence of these genes in our dataset (say 70%)
safe <- sapply(gos, function (x) {
  longi <- length(x)
  presence <- x %in% helaGenes
  percent <- sum(presence) / longi
  if (percent >= filtPercent) {
    TRUE
  } else {
    FALSE
  }
})

# Safe list of GOs
lgoshela <- names(gos[safe])


################# Pombe ##################
# Download GO data
# ftp://ftp.geneontology.org/pub/go/gene-associations/gene_association.pombase.gz
# fileUrl <- "ftp://ftp.geneontology.org/pub/go/gene-associations/gene_association.pombase.gz"
# download.file(fileUrl, 
#               destfile="./pombe/datos_cargar/gene_association.pombase.gz", 
#               method="curl")
# system("gunzip ./pombe/datos_cargar/gene_association.pombase.gz") 
pombeGO <- read.csv("./pombe/datos_cargar/gene_association.pombase", skip=44,
                      sep="\t", header=FALSE)
# Open HeLa dataset
pombe <- readRDS(file="pombe/datos_cargar/matriz_pombe.rds")
# pombe <- preprocessData(pombe)
pombe <- preprocessData2(pombe)
pombeGenes <- rownames(pombe)
# Try to tidy and reduce data
pos <- pombeGO[,2] %in% pombeGenes
pombeGO <- pombeGO[pos,]
lGOpombe <- pombeGO[,5]
lGOpombe <- as.character(lGOpombe)
lGOpombe <- sort(unique(lGOpombe))
# check number of ocurrences of each GO
safe <- sapply(lGOpombe, function (x) {
  pos <- pombeGO[,5] == x
  # continue if there are more than the chousen minimum
  if (sum(pos) >= mini) {
    genes <- unique(sort(as.character(pombeGO[pos,2])))
    if (length(genes) >= mini & length(genes <= maxi)) {
      TRUE
    } else {
      FALSE
    }
  } else {
    FALSE
  }
})
lGOpombe <- lGOpombe[safe]

# Now check for percentage 
safe <- sapply(lGOpombe, function (x) {
  pos <- pombeGO[,5] == x
  genes <- unique(sort(as.character(pombeGO[pos,2])))
  longi <- length(genes)
  presence <- genes %in% pombeGenes
  percent <- sum(presence) / longi
  if (percent >= filtPercent) {
    TRUE
  } else {
    FALSE
  }
})
lgosPombe <- lGOpombe[safe]


# Compare common presence of GO in the 3 species
pos1 <- lgosPombe %in% lgoshela
pos2 <- lgosPombe %in% lgosYeast
pos3 <- lgoshela %in% lgosYeast

# GO IDs in common between Pombe and HeLa
table(pos1)
# GO IDs in common between Pombe and saccharomyces
table(pos2)
# GO IDs in common between HeLa and saccharomyces
table(pos3)

# Get common GO terms to all 3 species
common <- lgosPombe[lgosPombe %in% lgoshela]
common <- lgosYeast[lgosYeast %in% common]

# Ahora toca sacar el termino al que se refieren en cada especie
# Get associated GO terms
GOtabl <- data.table(GO=common,
                     term=rep("NA", length(common)))
# We set GO as key
setkey(GOtabl, GO)
# start with HeLa
library("biomaRt")
GOtabl <- GoTermFinder(common, GOtabl, dataset="hsapiens_gene_ensembl", 
                       specie="term")

# NEED TO USE DIFFERENT METHODS, HOPE IS OK #
# Now we need a list of genes associated to each GOID and for each species
# hela
GOhela <- GoGeneFinder(common, dataset="hsapiens_gene_ensembl", 
                       id="hgnc_symbol")

# yeast (BioMart seems to lose much more genes)
gos <- org.Sc.sgdGO2ORF
gos <- as.list(gos)
GOyeast <- gos[common]

pos <- pombeGO[,5] %in% common
GOpombe <- pombeGO[pos, c(5,2)]




setwd(past_dir)
# PLOTEAR
# Get X edge names
nom_variables <- colnames(pombe)
nombre <- "plots_GO_comunes.pdf"
pdf(nombre)
# 2 rows and 3 columns per plot
par(mfrow=c(2,3))
# mai
# A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
for (com in common) {
  # yeast
  genesFilt_lev = unique(sort(as.character(GOyeast[[com]])))
  genesFilt_lev = as.character(genesFilt_lev[genesFilt_lev %in% rownames(lev)])
  # hela
  pos = GOhela[,1] == com
  genesFilt_hum = unique(sort(GOhela[pos,2]))
  genesFilt_hum = as.character(genesFilt_hum[genesFilt_hum %in% rownames(hela)])
  # pombe
  pos = GOpombe[,1] == com
  genesFilt_pomb = unique(sort(GOpombe[pos,2]))
  genesFilt_pomb = as.character(genesFilt_pomb[genesFilt_pomb %in% rownames(pombe)])
  
  
  # Calculamos medias y error estandar
  medias_lev = colMeans(lev[genesFilt_lev,], na.rm = TRUE)
  error_lev = apply(lev[genesFilt_lev,], 2, stderr)
  
  medias_hum = colMeans(hela[genesFilt_hum,], na.rm = TRUE)
  error_hum = apply(hela[genesFilt_hum,], 2, stderr)
  
  medias_pomb = colMeans(pombe[genesFilt_pomb,], na.rm = TRUE)
  error_pomb = apply(pombe[genesFilt_pomb,], 2, stderr)
  
  term = as.character(GOtabl[com, 2])
  # plot yeast
  titulo = paste(com, "_yeast\n", term, sep="")
  dim = 1:ncol(lev)
  plot(dim, medias_lev,
       ylim=range(c(0, 1)),
       pch=19, xlab="Variables", ylab="Mean +/- (SE)",
       main=titulo, xaxt = "n")
  lines(y = medias_lev[1:3], x = c(1,2,3))
  lines(y = medias_lev[4:6], x = c(4,5,6))
  axis(1, at=dim, labels=nom_variables, cex.axis=0.9)
  # hack: we draw arrows but with very special "arrowheads"
  arrows(dim, medias_lev-error_lev, dim, medias_lev+error_lev, length=0.05, angle=90, code=3)
  text(x=.85, y=0.005, labels=paste("maxN = ", length(genesFilt_lev), sep=""), pos = 4)
  
  # plot hela
  titulo = paste(com, "_HeLa\n", term, sep="")
  dim = 1:ncol(hela)
  plot(dim, medias_hum,
       ylim=range(c(0, 1)),
       pch=19, xlab="Variables", ylab="Mean +/- (SE)",
       main=titulo, xaxt = "n")
  lines(y = medias_hum[1:3], x = c(1,2,3))
  lines(y = medias_hum[4:6], x = c(4,5,6))
  axis(1, at=dim, labels=nom_variables, cex.axis=0.9)
  # hack: we draw arrows but with very special "arrowheads"
  arrows(dim, medias_hum-error_hum, dim, medias_hum+error_hum, length=0.05, angle=90, code=3)
  text(x=.85, y=0.005, labels=paste("maxN = ", length(genesFilt_hum), sep=""), pos = 4)
  
  # plot pombe
  titulo = paste(com, "_pombe\n", term, sep="")
  dim = 1:ncol(pombe)
  plot(dim, medias_pomb,
       ylim=range(c(0, 1)),
       pch=19, xlab="Variables", ylab="Mean +/- (SE)",
       main=titulo, xaxt = "n")
  lines(y = medias_pomb[1:3], x = c(1,2,3))
  lines(y = medias_pomb[4:6], x = c(4,5,6))
  axis(1, at=dim, labels=nom_variables, cex.axis=0.9)
  # hack: we draw arrows but with very special "arrowheads"
  arrows(dim, medias_pomb-error_pomb, dim, medias_pomb+error_pomb, length=0.05, angle=90, code=3)
  text(x=.85, y=0.005, labels=paste("maxN = ", length(genesFilt_pomb), sep=""), pos = 4)
}

dev.off()
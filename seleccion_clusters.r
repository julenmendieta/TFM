##seleccion_clusters.rmd
##2016-11-23 julenmendieta92@gmail.com
##1-Abrir los objetos de SOTA con 10, 15, 20 y 30 divisiones
##2-Tomar de cada uno los clusters indicados
##3-Generar un fichero con los GO asociados a cada cluster y otro con los plots de los mismos


############################################### Librerias  ########################################################
library("GO.db")
library("topGO")
library("GOstats")
library("biomaRt")
library("org.Sc.sgd.db")
library("org.Hs.eg.db")
library("GOSemSim")
library("tidyr")
library("clValid")
library("clusterProfiler")
library("openxlsx")

############################################### Funciones ######################################################
# Funcion para generar el plot
plotear <- function(DatosFilt, clu, clu_comp, ran1 = -0.1, ran2 = 1.1, go = TRUE, seguir = TRUE, valor_silueta, dir1=".", fin) {
  eje = 2
  for(i in 1: nrow(DatosFilt)) {
    if(i == 1) {
      plot(DatosFilt[1 ,], type ="l", col ="grey", ylim = range(c(ran1 ,ran2)),
           axes = FALSE , xaxt ="n", xlab ="", ylab ="")
      axis(1 , at =1:6 , labels = colnames(DatosFilt),
           cex.axis = eje)
    } else {
      par(new = TRUE)
      plot(DatosFilt[i ,] , type ="l", col ="grey", ylim = range(c(ran1 ,ran2)) ,
           axes = FALSE , xlab ="", ylab ="")
    }
  }
  par(new = TRUE)
  plot(apply(DatosFilt, 2, mean), type ="l", ylim = range(c( ran1 ,ran2)), xaxt = "n", 
       xlab ="Variables", ylab ="PromedioCluster", 
       main = paste("Cluster", "_", clu_comp, ".", as.character(clu), "_", nrow(DatosFilt), "_genes", sep =""),
       cex.axis = eje, sub=valor_silueta)
  lines(y = apply(DatosFilt[,1:3], 2, mean), x = c(1,2,3), col="blue", lwd=2)
  lines(y = apply(DatosFilt[,3:6], 2, mean), x = c(3,4,5,6), col="red", lwd=2)

  if ((min(DatosFilt) >= 0.35) & (max(DatosFilt) <= 0.65) & (seguir == TRUE)) {
    plotear(DatosFilt, clu, clu_comp, ran1 = 0.35, ran2 = 0.65, go = FALSE, seguir = FALSE)
  }
  if (go == TRUE) {
    salida = enriquecimientoGO(rownames(DatosFilt), fin)
    # Reducimos los decimales que aparecen para los GO en SOTA
    # options("scipen"=-100, "digits"=3)
#     for (nom in names(salida)) {
#       # Reducimos la longitud del nombre del GO si es demasiado largo
#       if (length(salida[[nom]]) > 1) {
#         for (c in 1:nrow(salida[[nom]])){
#           if (nchar(salida[[nom]][c, "Term"]) > 74) {
#             salida[[nom]][c, "Term"] = paste(substr(salida[[nom]][c, "Term"], 1, 74), "...", sep = "")
#           }
#         }
#       }
#       print(paste("#############################", nom, "######################################"))
#       print(as.data.frame(salida[nom]), row.names = FALSE)
#     }
    # Volvemos a las condiciones por defecto, que si no todos los numeros se muestran asi
    # options("scipen"= NULL, "digits"=7)
    # Guardamos en un fichero las GO
    write.xlsx(salida, file = paste(dir1, "/salida_GO/", "Cluster", "_", clu_comp, ".", as.character(clu), 
                                      "_", ... = nrow(DatosFilt), "_genes.xlsx", sep=""))
  }
}


# Enriquecimineto GO
enriquecimientoGO <- function(listaGenes, fin) {
  ontologies = c("BP", "CC")
  if (fin == "_levadura") {
      orfId <- org.Sc.sgdGO
      anotacion = "org.Sc.sgd"
      especie = "yeast"
      
      universe <- mappedkeys(orfId)
      catego1 <- as.character(listaGenes)
  } else if (fin == "_HeLa") {
      orfId <- org.Hs.egGO
      anotacion = "org.Hs.eg"
      especie = "human"
      
      # Cambiamos los IDs a entrez
      xx1 <- as.list(org.Hs.egALIAS2EG)
      listaGenes1 = xx1[listaGenes]
      names(listaGenes1) = NULL
      catego1 = unlist(listaGenes1)
      catego1 = unique(catego1)
      
      # Guardamos una lista con todos los genes que tienen GO asociados
      orfId <- org.Hs.egGO
      conGO <- mappedkeys(orfId)
      # Guardamos una lista con todos los genes
      universe <- xx1
      names(universe) = NULL
      universe = unlist(universe)
      universe = unique(universe)
      
  }
  
  
  salida = list()
  for(ont in ontologies) {
    params <- new("GOHyperGParams" , geneIds = catego1,
                universeGeneIds = universe , ontology = ont , pvalueCutoff =0.001,
                conditional =F , testDirection = "over" , annotation = anotacion)
    hgOver <- hyperGTest(params)
    
    
    ## Ajustamos el p-valor con FDR
    # Get the p-values of the test
    hgOver.pv <- pvalues(hgOver)
    # Adjust p-values for multiple test (FDR)
    hgOver.pv.fdr <- p.adjust(hgOver.pv,'fdr')
    # Guardamos los valores del FDR en el objeto de resultados
    result <- summary(hgOver)
    # Si no hay resultados, esta parte sobra
    if(dim(result)[1] != 0) {
      longi = dim(result)[1]
      FDR = hgOver.pv.fdr[1 : longi]
      names(FDR) = NULL
      result["FDR"] <- FDR
      if(ont == "BP") {
        resultsalida <- result[, c("Term", "FDR", "GOBPID")]
      } else {
        resultsalida <- result[, c("Term", "FDR", "GOCCID")]
      }
      # select the GO terms with adjusted p-value less than the cut off
      
      resultsalida = simplificar(resultsalida, ont, cutoff=0.7, by="FDR", 
                                 select_fun=min, measure="Rel", organismo = especie)
      

      resultsalida <- resultsalida[resultsalida[, "FDR"] < 0.001, ]
      if (nrow(resultsalida) == 0) {
        salida[[ont]] = "No summarized results with FDR < 0.001"
      } else {
        salida[[ont]] = resultsalida
      }
    } else {
      salida[[ont]] = "There weren't adjusted p values < 0.001 before GO summary"
    }
  }
  
  return(salida)
}


# Funcion para eliminar GO solapantes por similitud lexica (0.7)
simplificar <- function(x, ont, cutoff=0.7, by="FDR", select_fun=min, measure="Rel", organismo = "otro") {
      ## to satisfy codetools for calling gather
      # Con esto indicamos la similitud entre cada GO
      go1 <- go2 <- similarity <- NULL
      
      ID = paste("GO", ont, "ID", sep="")
      rownames(x) = x[,ID]
      
      
      sim <- mgoSim(x[,ID], x[,ID],
                    ont=ont,
                    organism=organismo,
                    measure=measure,
                    combine=NULL)
      
      sim.df <- as.data.frame(sim)
      sim.df$go1 <- row.names(sim.df)
      sim.df <- gather(sim.df, go2, similarity, -go1)
      
      sim.df <- sim.df[!is.na(sim.df$similarity),]
      
      ## feature 'by' is attached to 'go1'
      sim.df <- merge(sim.df, x[, c(ID, by)], by.x="go1", by.y=ID)
      sim.df$go2 <- as.character(sim.df$go2)
      
      ID2 <- x[,ID]
      
      GO_to_remove <- character()
      for (i in seq_along(ID2)) {
          ii <- which(sim.df$go2 == ID2[i] & sim.df$similarity > cutoff)
          ## if length(ii) == 1, then go1 == go2
          if (length(ii) < 2) 
              next
          
          sim_subset <- sim.df[ii,]
          
          jj <- which(sim_subset[, by] == select_fun(sim_subset[, by]))
          
          ## sim.df <- sim.df[-ii[-jj]]
          GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
      }
      x <- x[!(x[,ID] %in% GO_to_remove),]
      return(x)
}

# Funcion para generar silueta y lanzar las funciones para generar plot y guardar GOs
selecClusters <- function(div_clusters1, dir1, grupos1, go=TRUE, fin) {
    for ( i in 1:length(div_clusters1)) {
        s.sota = readRDS(paste(dir1, "/salida",fin, "_", div_clusters1[i], "div", "/0", sep=""))
        Matriz_filtr = s.sota$data
        # Analizamos la robustez de los clustes generados con silhouette
        distance = dist(Matriz_filtr, method = "euclidean")
        clusters_int = as.integer(s.sota$clust)
        names(clusters_int) = rownames(Matriz_filtr)
        silueta = silhouette(x=clusters_int, dist=distance)
        resumen = summary(silueta)
        silueta_clust = resumen$clus.avg.widths
        silueta_clust = round(silueta_clust, 2)
        for (clu in grupos1[,i]) {
            if (!is.na(clu)) {
                pos = s.sota$clust == clu
                DatosFilt = Matriz_filtr[pos,]
                valor_silueta = silueta_clust[as.numeric(clu)]
                clu_comp = paste(div_clusters1[i], "div", "_0",sep="")
                plotear(DatosFilt, clu, clu_comp, go=go, valor_silueta=valor_silueta, dir1=dir1, fin=fin)
            }
        }

    }
}


# La idea es generar un excel para cada division, y con los clusters separados por hojas. Todo esto a partir de la salida
# del paso anterior. Estamos agrupando todos los excel que se han generado
juntarFicheros <- function (dir1, div_clusters, grupos) {
    dir2 = paste(dir1, "/salida_GO/", sep="")
    ficheros = list.files(dir2)
    tabla_GO = list()
    # colnames(tabla_GO) = paste(div_clusters, "div", sep="")
    for ( i in div_clusters) {
        tabla_temp = list()
        pos_fich = grep(paste("_", i, "div", sep=""), ficheros)
        division = ficheros[pos_fich]
        # Vamos abriendo cada fichero
        for (gru in 1:length(grupos[,as.character(i)])) {
            # Si el cluster es NA pasamos
            if (!is.na(division[gru])) {
                # Abrimos la primera hoja y modificamos los nombres de columna para la ontologia a la que se refiere
                fich_temp = read.xlsx(paste(dir2, division[gru], sep=""), sheet = "BP")
                fich_temp2 = fich_temp
                columnas1 = ncol(fich_temp)
                nom_columnas1 = colnames(fich_temp)
                colnames(fich_temp2) = paste(nom_columnas1, "_BP", sep="")
                # Le anyadimos unas columnas extra para separar de la siguiente ontologia
                # Si no hay GOs da error
                if (columnas1 != 1) {
                    fich_temp2[,(columnas1+1):(columnas1+3)] = ""
                }
                # Abrimos la segunda hoja y los juntamos
                fich_temp = read.xlsx(paste(dir2, division[gru], sep=""), sheet = "CC")
                # Para ello, primero modificamos los nombres de columnas
                nom_columnas2 = colnames(fich_temp)
                columnas2 = ncol(fich_temp)
                colnames(fich_temp) = paste(nom_columnas2, "_CC", sep="")
                # Luego igualamos el numero de filas
                lin1 = nrow(fich_temp2)
                lin2 = nrow(fich_temp)
                if (lin1 > lin2) {
                    resto = lin1 - lin2
                    fich_temp[(lin2 + 1):(lin2 + resto),] = ""
                } else if (lin1 < lin2) {
                    resto = lin2 - lin1
                    fich_temp2[(lin1 + 1):(lin1 + resto),] = ""
                }
                fich_temp = cbind(fich_temp2, fich_temp)
                # Si en uno de las subontologias habia GO y en la otra no, solo se habra guardado una
                if (length(nom_columnas1) == 1 & length(nom_columnas2) != 1) {
                    fich_temp[,(columnas2+1):(columnas2+3)] = ""
                    colnames(fich_temp) = c(colnames(fich_temp)[1:5],
                                            paste("No.summarized.results.with.FDR.<.0.001", "_BP", sep="" ))
                } else if (length(nom_columnas1) != 1 & length(nom_columnas2) == 1) {
                    colnames(fich_temp) = c(colnames(fich_temp)[1:5],
                                            paste("No.summarized.results.with.FDR.<.0.001", "_CC", sep="" ))
                    
                }
                
                # Anyadimos esta tabla a la tabla asociada a esta division
                clu = as.character(grupos[,as.character(i)][gru])
                tabla_temp[[clu]] = fich_temp
            }
        }
        # Guardamos esta lista en otra
        div_clu = paste(i, "div", sep = "")
        tabla_GO[[div_clu]] = tabla_temp
    }
    
    # Generamos un excel por cada una de las divisiones
    dir.create(file.path(paste(dir1, "/salida_GO_juntada", sep="")), showWarnings = FALSE)
    for (divi in names(tabla_GO)) {
        temporal = tabla_GO[[divi]]
        write.xlsx(temporal, file = paste(dir1, "/salida_GO_juntada/", "Terminos_GO_", divi, ".xlsx", sep=""))
    }
}


################################################# Programa  #####################################################
# # levadura
# diez = c(2,4,6,8,10,11)
# quince = c(8,13,14,11,9,10)
# veinte = c(14,15,NA,NA,NA,NA)
# treinta = c(30,31,16,17,NA,NA)

# dir1 = "/media/julen/TOSHIBA_EXT/datos/paper/levadura"
# fin = "_levadura"
# fich_salida = "Selección_clusters_levadura.pdf"

# HeLa
diez = c(5,8,9,10,NA,NA,NA,NA,NA,NA,NA,NA)
quince = c(13,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
veinte = c(6,12,13,16,17,20,21,NA,NA,NA,NA,NA)
treinta = c(12,13,14,15,20,21,22,23,24,25,26,27)

dir1 = "/media/julen/TOSHIBA_EXT/datos/paper/HeLa"
fin = "_HeLa"
fich_salida = "Selección_clusters_HeLa.pdf"


### Parte fija  ###

setwd(dir1)

# Seleccion de los clusters
div_clusters = c(10,15,20,30)
grupos <- matrix(c(diez, quince, veinte, treinta),ncol=4)
colnames(grupos) = div_clusters

# Creamos la carpeta de salida de los GO
dir.create(file.path(paste(dir1, "/salida_GO", sep="")), showWarnings = FALSE)

# Lanzamos el programa
pdf(fich_salida)
par(mfrow=c(2,2))
selecClusters(div_clusters, dir1, grupos, go=FALSE, fin)
dev.off()


# Ahora tenemos que juntar todos los excel generados
juntarFicheros(dir1, div_clusters, grupos)



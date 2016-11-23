##cluster_heatmap.r
##2016-11-23 julenmendieta92@gmail.com
##1-Abrir matriz de datos de expresion genica, imputar datos faltantes y convertir a ranking
##2-Generar HeatMaps utilizando distancia euclidea y de manhattan y el metodo UPGMA
##3-Comprobar robustez del clustering mediante silhouette

dir = "/media/julen/TOSHIBA_EXT/datos/paper/levadura"
setwd(dir)
matriz = readRDS(file="datos_cargar/matriz_lev.rds")
fin="_levadura"

##################### FUNCIONES ##################################
# Vamos a plantear una funcion con silhouette que nos diga en que punto los clusters estan mejor hechos
# (puede ser que encontremos la mejor situacion de un clustering que no vale para nada)
buscarSilhouette <- function(heat1, dir, distan, distance) {
  # Obtenemos la h del dendrograma para ver hasta donde separar clusters
  height = max(cophenetic(heat1$rowDendrogram))
  # Este valor esta escogido despues de ver donde empiezan a aglutinarse las divisiones en los dos arboles
  desde = height/2.675
  mejor = c(0,0)
  maxi = c(0,0)
  for (prueba in seq(from=desde,to=(height-0.1), length.out = 12)) {
    # El valor minimo razonable es de 0.51, dicho esto, me voy a quedar con el punto en el que mas clusters superen
    # este valor. Ojo, porque puede ser que perdamos por el camino un cluster muy bueno
    
    # Primero dividimos el cluster, usamos silhoutte y miramos cuantos clusters tienen un valor aceptable
    clusters = cutree(heat1$rowDendrogram, h=prueba)
    table(clusters)
    silueta = silhouette(x=clusters, dist=distance)
    resumen = summary(silueta)
    resumen = resumen$clus.avg.widths
    filtro_ok = table(resumen > 0.5)
    
    # Si tenemos mas clusters pasando el filtro que antes, guardamos este nivel de corte
    if (!is.na(filtro_ok[2]) & filtro_ok[2] > mejor[2]) {
      mejor[2] = filtro_ok[2]
      mejor[1] = prueba
    }
    # Tambien vamos a buscar el cluster con mayor valor de silhouette, para guardarlo
    if (max(resumen) > maxi[2]) {
      maxi = c(prueba, max(resumen))
    }
  }
  # Ahora que tenemos esta seleccion, pasamos la informacion de esos clusters con los genes a un fichero
  # Primero guardamos la lista del punto en el que mas clusters cumplian con el requisito minimo
  clusters = cutree(heat1$rowDendrogram, h=mejor[1])
  silueta = silhouette(x=clusters, dist=distance)
  resumen = summary(silueta)
  resumen = resumen$clus.avg.widths
  # Guardamos solo los clusters con disimilaridad media superior a 0.5
  resumen = resumen[resumen > 0.5]
  clusters_top = data.frame()
  for (clu1 in names(resumen)) {
    genes = names(clusters[clusters == clu1])
    clusters_top[1:length(genes), paste(clu1, round(resumen[[clu1]], digits = 2), sep = "_")] = genes
  }
  
  # Ahora pasamos al cluster con mejor resultado de entre todos
  clusters = cutree(heat1$rowDendrogram, h=maxi[1])
  silueta = silhouette(x=clusters, dist=distance)
  resumen = summary(silueta)
  resumen = resumen$clus.avg.widths
  resumen = resumen[resumen == maxi[2]]
  cluster_maxi = data.frame()
  genes = names(clusters[clusters == names(resumen)])
  cluster_maxi[1:length(genes), paste(names(resumen), round(resumen[[names(resumen)]], digits = 2), sep = "_")] = genes
  
  # Guardamos
  nombre = paste(dir,"/heatMap/", distan,"/MejoresClusters_silhouette_h=", round(mejor[1],digits=2), sep="")
  write.table(clusters_top, file=nombre, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  nombre = paste(dir,"/heatMap/", distan,"/MejorClustersGeneral_silhouette_h=", round(maxi[1], digits=2), sep="")
  write.table(cluster_maxi, file=nombre, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  # table(clusters[clusters %in% names(resumen[resumen > 0.5])])
}

###################################################### LIBRERIAS #################################################
library("impute")
library("cluster")
library(gplots)
library(RColorBrewer)
library(dendextend)
library(colorspace) # get nice colors

###################################################### CODIGO ###################################################



matriz = matriz[,c("TR", "RS", "RA", "TLRi", "PS", "PA")]

# Guardamos solo los genes de los que tengamos info en al menos 4 variables, y de distintas moleculas
# Para SOTA utilizamos TLRi
matriz1 = matriz[,c("TR","RS","RA")]
posiciones= apply(!is.na(matriz1), 1, sum)
pos = posiciones >= 2
print("Genes con 2 de 3 variables en RNA")
table(pos)
matriz1 = matriz1[pos,]

temporal = impute.knn(as.matrix(matriz1) ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
matriz1 = as.data.frame(temporal$data)

matriz2 = matriz[,c("TLRi","PS","PA")]
posiciones= apply(!is.na(matriz2), 1, sum)
pos = posiciones >= 2
print("Genes con 2 de 3 variables en RNA")
table(pos)
matriz2 = matriz2[pos,]
# Generamos valores faltantes con el knn
temporal = impute.knn(as.matrix(matriz2) ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
matriz2 = as.data.frame(temporal$data)

# Juntamos las dos matrices
pos = rownames(matriz1) %in% rownames(matriz2)
table(pos)
Matriz_6 = data.frame(row.names = rownames(matriz1)[pos])
Matriz_6[,c("TR","RS","RA")] = matriz1[rownames(Matriz_6),c("TR","RS","RA")]
Matriz_6[,c("TLRi","PS","PA")] = matriz2[rownames(Matriz_6),c("TLRi","PS","PA")]

# Pasamos a matriz asegurandonos de que si alguna columna esta como factor, no se guarde su indice
indx <- sapply(Matriz_6, is.factor)
if (sum(indx) > 0) {
  Matriz_6[indx] <- lapply(Matriz_6[indx], function(x) as.numeric(as.character(x)))
}
Matriz_6 = as.matrix(Matriz_6)

# Generamos ranking
for(col in 1:ncol(Matriz_6)) {
  Matriz_6[, col] = rank(Matriz_6[, col], na.last = "keep")
}
# Lo pasamos a escala de 0 a 1
Matriz_6 = t(t(Matriz_6)/apply(Matriz_6,2,max, na.rm = TRUE))


# Computes agglomerative hierarchical clustering of the dataset.
# Se permiten NAs
# arbol = agnes(Matriz_6, diss = FALSE, metric = "euclidean",
      # stand = FALSE, method = "average", keep.data = TRUE, trace.lev = 0)




################ Heat map #######################################


mat_data = Matriz_6

print("Para el generar los clusters, la matriz de distancias se calcula mediante el metodo euclideo,
      y el clustering mediante UPGMA")
suppressWarnings(dir.create(file.path (dir, "heatMap")))
# Para el clustering vamos a plantearlo con distancia eclidea y de manhatan
for (distan in c("euclidean", "manhattan")) {
  suppressWarnings(dir.create(file.path (dir, "heatMap", distan)))
  # method	
  # the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", 
  # "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
  distance = dist(mat_data, method = distan)
  clustering = hclust(distance, method = "average")
  
  # creates a own color palette from red to green
  my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
  
  # (optional) defines the color breaks manually for a "skewed" color transition
  col_breaks = c(seq(0,0.33,length=100),  # for red
                 seq(0.34,0.66,length=100),           # for yellow
                 seq(0.67,1,length=100))             # for green
  
  # # creates a 5 x 5 inch image
  # png("./heatMap/heatmaps_in_r.png",    # create PNG for the heat map        
  #     width = 5*600,        # 5 x 300 pixels
  #     height = 5*600,
  #     res = 600,            # 300 pixels per inch
  #     pointsize = 8)        # smaller font size
  
  # Primero generamos el heatmap a lo grande
  jpeg(filename = paste(dir,"/heatMap/", distan, "/heatMap_lev.jpeg", sep=""),
       width = 4*600, height = 4*600, units = "px", pointsize = 15,
       quality = 100,
       bg = "white")
  
  heat = heatmap.2(mat_data,
            main = "HeatMap global", # heat map title
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(12,9),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            breaks=col_breaks,    # enable color transition at specified limits
            dendrogram="row",     # only draw a row dendrogram
            Colv=FALSE,           # turn off column clustering
            Rowv = as.dendrogram(clustering)) # Con esto indicamos metodo de clustering en las filas
  
  
  dev.off()               # close the PNG device
  
  ## Vamos a enseñar lo que hemos dividido
  # Primero generamos un objeto en el que dividimos los genes en div grupos
  if (distan == "manhattan") {
    div = 5
  } else {
    div=6
  }
  clusters = cutree(heat$rowDendrogram,k=div)
  # Luego creamos una matriz con los genes y el cluster al que se han asociado
  mat_genesClust = matrix(append(names(clusters), rep(NA,length(clusters))), ncol=2)
  rownames(mat_genesClust) = mat_genesClust[,1]
  # Asociamos los nombres de los genes a cada uno de esos clusters
  for (c in 1:div) {
    clust = clusters[clusters == c]
    clust = names(clust)
    mat_genesClust[clust,2] = c
  }
  # Preparamos el dendrograma
  dend <- rotate(heat$rowDendrogram, 1:dim(mat_data)[1])
  # Coloreamos las ramas en funcion de las divisiones que hemos planteado
  dend1 <- color_branches(dend, k = div)
  # Asociamos los nombres de los genes a cada uno de esos clusters para el plot
  labels_colors(dend1) <-
    rainbow_hcl(div)[sort_levels_values(
      as.numeric(mat_genesClust[,2])[order.dendrogram(dend)]
    )]
  # Reducimos el tamaño de las labels 
  dend1 <- set(dend1, "labels_cex", 0.5)
  pdf(file = paste(dir,"/heatMap/", distan, "/Division_cluster_inicial.pdf", sep=""))
  plot(dend1, main="Division de los clusters para HeatMap individuales", 
       sub = paste("nº de genes por clusters ordenados=", paste(as.character(table(clusters)), collapse = ", "), sep=""))
  legend("topleft", legend = unique(mat_genesClust[order.dendrogram(dend),2]), fill = rainbow_hcl(div))
  dev.off()               # close the JPEG device
  
  
  # Miramos la robusted de esas divisiones con silhouette
  buscarSilhouette(heat, dir, distan, distance)
  
  ## Ahora toca generar los heatMaps individuales
  table(clusters)
  nombre = "./heatMap/heatMap_lev_"
  for(n in 1:div) {
    # Por cada cluster tomamos los ID de los genes
    clust = clusters[clusters == n]
    clust = names(clust)
    mat_data2 = Matriz_6[clust,]
    # Preparamos el metodo de clustering
    distance = dist(mat_data2, method = distan)
    clustering = hclust(distance, method = "average")
    # Preparamos la longitud del PDF en funcion del numero de genes
    longi = dim(mat_data2)[1]
    if(longi <= 10) {
      alt=10
      anch=6
      l1=0.3
      l2=0.7
      mar2= 12
    } else if(longi < 100 && longi > 10){ 
      alt=25
      anch=5
      l1=0.03
      l2=1
      mar2=5
    } else if(longi >= 100 && longi < 1100){
      alt=50
      anch=10
      l1=0.1
      l2=1
      mar2=5
    } else {
      alt=70
      anch=10
      l1=0.03
      l2=1
      mar2=5
    }
    nombre_temp = paste(dir,"/heatMap/", distan, "/heatMap_lev", n, ".pdf", sep="")
    # Generamos otro heat map con esos genes
    pdf(file=nombre_temp, height=alt, width=anch)
    heat_temp = heatmap.2(mat_data2, trace='none', dendrogram='row', Colv=F, scale='none', 
              col=my_palette, symbreak=T,
              margins=c(5,mar2), keysize=1,
              lwid=c(1,0.05,1), lhei=c(l1,l2), lmat=rbind(c(5,0,4),c(3,1,2)),
              RowSideColors=as.character(clusters), breaks = seq(0, 1, length.out = 300),
              Rowv = as.dendrogram(clustering))
    
    dev.off()               # close the PDF device
    
    # Ahora necesitamos una lista con los genes en orden para cada cluster generado
    genes = rownames(mat_data2[order.dendrogram(heat_temp$rowDendrogram),])
    # Como la lista viene al reves le damos la vuelta
    genes = rev(genes)
    write.table(genes, file=paste(dir,"/heatMap/", distan,"/Genes_cluster_", n, sep=""), 
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}


# Originales
# heatmap.2(mat_data2, trace='none', dendrogram='row', Colv=F, scale='none', 
#           hclust=hclustfunc, distfun=distfunc, col=greenred(256), symbreak=T,
#           margins=c(10,20), keysize=0.5,
#           lwid=c(1,0.05,1), lhei=c(0.03,1), lmat=rbind(c(5,0,4),c(3,1,2)),
#           RowSideColors=as.character(clusters))

# heatmap.2(mat_data2,
#           main = "Correlation", # heat map title
#           density.info="none",  # turns off density plot inside color legend
#           trace="none",         # turns off trace lines inside the heat map
#           margins =c(12,9),     # widens margins around plot
#           col=my_palette,       # use on color palette defined earlier
#           breaks=col_breaks,    # enable color transition at specified limits
#           dendrogram="row",     # only draw a row dendrogram
#           Colv="NA",
#           cexCol = 5)            # turn off column clustering








############## NOTAS
# El valor "desde" en buscarSilhouette() se tiene que cambiar en funcion de los arboles que den cada metodo
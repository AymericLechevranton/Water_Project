# Projet Classification - M1 Ingénierie Statistique - Lechevranton & Vasse

# Introduction : Ce script s'organise en 2 parties.
# Partie 1 : la mise en oeuvre et la comparaison de méthodes de classification non supervisée.
# Partie 2 : la mise en oeuvre de méthodes de classification supervisée pour prévoir la nature de l'eau.



######################## Partie 0 - Travail préliminaire ######################## 

    # Importer les librairies nécessaires
library(cluster)
library(readr)
library(factoextra)
library(ade4)
library(FactoMineR)
library(MASS)
library(dplyr)
library(KernSmooth)
library(caret)
#source("https://bioconductor.org/biocLite.R") ; biocLite("ComplexHeatmap")
library(ComplexHeatmap)
library(clustMixType)
library(NbClust)
library(clValid)
library(fpc)
library(devtools)
#install_github("larmarange/JLutils")
library(JLutils)
library(RColorBrewer)
library(formattable)
library(rpart)
library(dbscan)
library(bootstrap)
library(clusterCrit)
library(psych)
library(sortinghat)
library(class)
library(fmsb)
library(mvnormtest)
library(heplots)




    # Importer les données nécessaires
eaux <- data.frame(read_delim("C:/Users/Aymeric/Documents/#M1 - Université/M1-S1_-_Analyse de données/Projet/Eaux2018.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE))
eaux2018 <- eaux
head(eaux2018)

#Afficher le tableau du jeu données utilisé
datatable<-formattable(data.frame(na.omit(eaux2018)[,3:11]))
datatable

    # Modification de la variable qualitative 'Nature' et 'Pays'
eaux2018$Nature <- ifelse(eaux2018$Nature=='gaz', 1, 0)
eaux2018$Pays <- ifelse(eaux2018$Pays=='France', 1, 0)
head(eaux2018)

    # On ne choisit que les données "complètes" (ie : sans NA) et nous les centrons réduisons.
data <- scale(na.omit(eaux2018)[,3:11])
head(data)

    #dim du jeu de données
dim(data)
summary(data)
boxplot(data)
#62 individus et 9 variables quantitatives



      # Statistique de Hopkins pour connaître la tendance de regroupement des données.
get_clust_tendency(data,2)$hopkins_stat #stat de Hopkins = 0.165 (<0.5)
#donc le jeu de données est intéressant à "clusturiser".
get_clust_tendency(data,2)$plot+labs(title ="Vérification graphique - VAT")
# Nous obervons bien qu'il existe bien une tendance au regroupement dans le jeu de données.



######################## Partie 1 - Classification non supervisée ######################## 

# Dans cette partie, nous allons utiliser différentes méthodes de classification non supervisée :
# - Classification par partition et classification floue (centres mobiles (Forgy), kmeans (MacQueen), Noyau/Nuées dynamiques (Hartigan-Wong), k-medoïdes (PAM et CLARA), FuzzK (FANNY)).
# - Classification Hiérarchique (ascendante (AGNES) et descendante (DIANA)).
# - Classification par voisinage dense (kernel).
# - La classification mixte ici n'a pas d'intérêt car nous avons peu de données.
# - Double CAH.

    


    #################### Méthodes de classification ####################


    # Classification kmeans

# Nombre optimal de clusters pour k-means : choix arbitraire.

# average silhouette ("silhouette"), within cluster sums of squares ("wss"), gap statistics ("gap_stat").
fviz_nbclust(data, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
# k=2
fviz_nbclust(data, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method") 
# k=3
fviz_nbclust(data, kmeans, method = "gap_stat")+
  labs(subtitle = "Gap statistic method")
# k=2

# Par la suite, nous prendrons arbitrairement k=2.

#'Hartigan-Wong' (nuées dynamiques par défaut), 'Lloyd', 'Forgy' (centres mobiles) et 'MacQueen' (kmeans).
km.res <- kmeans(data, 2, algorithm='Hartigan-Wong')

#Représentation fviz
fviz_cluster(km.res, data = data,
             palette = c("#00AFBB","#2F8F00",
                         "#2E9FDF","#259FEF",
                         "#006600"),
             ggtheme = theme_minimal(),
             main = "Graphique de partitionnement kmeans"
)

#Autres représentations
s.class(data,as.factor(km.res$cluster), col = brewer.pal(3, "Set1"), sub = "Axes 1 et 2")
clusplot(data,as.factor(km.res$cluster))



    # Classification k-médoïdes (PAM - Partition Around Medoids)

# Nombre optimal de clusters pour k-medoids/pam : choix arbitraire.

# average silhouette ("silhouette"), within cluster sums of squares ("wss"), gap statistics ("gap_stat").
fviz_nbclust(data, pam, method = "silhouette")+
  labs(subtitle = "Silhouette method") 
# k=3
fviz_nbclust(data, pam, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method") 
# k=3
fviz_nbclust(data, pam, method = "gap_stat")+
  labs(subtitle = "Gap statistic method")
# k=10

# Par la suite, nous prendrons arbitrairement k=3.

pam.res <- pam(data, 3)

#Représentation fviz
fviz_cluster(pam.res, data = data,
             palette = c("#00AFBB","#2F8F00",
                         "#2E9FDF","#259FEF",
                         "#006600"),
             ggtheme = theme_minimal(),
             main = "Graphique de partitionnement PAM"
)

#Autre représentation
s.class(data,as.factor(pam.res$cluster), col = brewer.pal(4, "Set1"), sub = "Axes 1 et 2")
plot(pam.res, which=1)
plot(pam.res, which=2)

#Nous obervons que la silhouette moyenne est de 0.52 (62), 
#avec 0.59 (44) pour la partition 1, 0.21 (14) pour la partition 2 et 0.78 (4) pour la partition 3.



    # Classification k-médoïdes (CLARA - Clustering LARge Application)

# L'algorithme PAM demande beaucoup de calculs. L'algorithme CLARA est préféré à PAM dans le cas de données importantes.

# Nombre optimal de clusters pour k-medoids/pam : choix arbitraire.

# average silhouette ("silhouette"), within cluster sums of squares ("wss"), gap statistics ("gap_stat").
fviz_nbclust(data, clara, method = "silhouette")+
  labs(subtitle = "Silhouette method") 
# k=3
fviz_nbclust(data, clara, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method") 
# k=3
fviz_nbclust(data, clara, method = "gap_stat")+
  labs(subtitle = "Gap statistic method")
# k=10

# Par la suite, nous prendrons arbitrairement k=3.

clara.res <- clara(data, 3)

#Représentation fviz
fviz_cluster(clara.res, data = data,
             palette = c("#00AFBB","#2F8F00",
                         "#2E9FDF","#259FEF",
                         "#006600"),
             ggtheme = theme_minimal(),
             main = "Graphique de partitionnement CLARA"
)

#Autre représentation
s.class(data,as.factor(clara.res$cluster), col = brewer.pal(4, "Set1"), sub = "Axes 1 et 2")
plot(clara.res,which=1)
plot(clara.res,which=2)

#Nous obervons que la silhouette moyenne est de 0.53 (46), 
#avec 0.63 (29) pour la partition 1, 0.21 (13) pour la partition 2 et 0.79 (4) pour la partition 3.


    

    # Classification floue (FANNY)

# Nombre optimal de clusters pour FANNY : choix arbitraire.

# average silhouette ("silhouette"), within cluster sums of squares ("wss"), gap statistics ("gap_stat").
fviz_nbclust(data, fanny, method = "silhouette")+
  labs(subtitle = "Silhouette method") 
# k=2
fviz_nbclust(data, fanny, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2)+
  labs(subtitle = "Elbow method") 
# k=2
fviz_nbclust(data, fanny, method = "gap_stat")+
  labs(subtitle = "Gap statistic method")
# k=6
# Par la suite, nous prendrons arbitrairement k=2.

fanny.res <- fanny(data, 2)
fanny.res$coeff
# Coefficient de Dunn de 0.548 (de 0 à 1 : de "completely fuzzy" à "hard cluster").

#Représentation fviz
fviz_cluster(fanny.res, data = data,
             palette = c("#00AFBB","#2F8F00",
                         "#2E9FDF","#259FEF",
                         "#006600"),
             ggtheme = theme_minimal(),
             main = "Graphique de partitionnement FANNY"
)

#Autre représentation
s.class(data,as.factor(fanny.res$cluster), col = brewer.pal(2, "Set1"), sub = "Axes 1 et 2")
plot(fanny.res, which=1)
plot(fanny.res, which=2)

#Nous obervons que la silhouette moyenne est de 0.42 (62), 
#avec 0.65 (40) pour la partition 1 et -0.009 (22) pour la partition 2.



    # Classification ascendante hiérarchique (AGNES - AGglomerative NESting)

  # Choix de la méthode pour AGNES

data <- scale(na.omit(eaux2018)[,3:11])
head(data)

results <- function(j){
corcoph <- data.frame(matrix(0, nrow=4, ncol=6))
colnames(corcoph) <- c('average', 'single', 'complete', 'ward.D','ward.D2','mcquitty')
rownames(corcoph) <- c('Corrélation cophénétique', 'Silhouette', 'Indice de Dunn', 'Connectivité')
for (i in c('average', 'single', 'complete', 'ward.D','ward.D2','mcquitty')){
  hc <- eclust(data, FUNcluster="hclust", hc_method=i, k=j)
  corcoph[1,i] <- round(cor(dist(data),cophenetic(hc)), 3)
  corcoph[2,i] <- hc$silinfo$avg.width
  corcoph[3,i] <- cluster.stats(dist(data), hc$cluster)$dunn
  corcoph[4,i] <- connectivity(dist(data), hc$cluster)
}
corcoph <- round(corcoph,3)
# Matrice des points obtenus au classement des méthodes
agnes_choice <- rank(corcoph[1,],ties.method="max")+rank(corcoph[2,],ties.method="max")+rank(corcoph[3,],ties.method="max")+(6-rank(corcoph[4,],ties.method="max"))
return(data.frame(t(agnes_choice)))
}

res.model <- t(sapply(2:9,results))
rownames(res.model) <- c('k=2','3','4','5','6','7','8','9')
res.model <- formattable(data.frame(res.model))
res.model
# single et complete sont les meilleurs méthodes :
# nous choisissons complete qui nous permet d'avoir un bon modèle pour un petit k, ce qui nous permettra d'interpréter plus facilement.

  # Nombre optimal de clusters pour AGNES : choix arbitraire.

data <- scale(na.omit(eaux2018)[,3:11])
agnes.res <- hclust(dist(data), method='complete')

#Sauts d'inertie
inertie <- sort(agnes.res$height, decreasing = TRUE)
plot(inertie[1:20], type = "s", xlab = "Nombre de classes", ylab = "Inertie", main= "Graphique des sauts d'inertie")
points(c(2,3,6), inertie[c(2,3,6)], col = c("green3","red3","blue3"), cex = 2, lwd = 3)
#k={2,3,6} se distinguent des autres

#Calcul de l'inertie relative
best.cutree(agnes.res, graph = TRUE, xlab = "Nombre de classes", ylab = "Inertie relative")
#La classification k=3 est devant toutes les autres.

# average silhouette ("silhouette"), within cluster sums of squares ("wss"), gap statistics ("gap_stat").
fviz_nbclust(data, hcut, hc_func = "hclust", hc_method = "complete", method = "silhouette")+
  labs(subtitle = "Silhouette method") 
# k=4
fviz_nbclust(data, hcut, hc_func = "hclust", hc_method = "complete", method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method") 
# k=4
fviz_nbclust(data, hcut, hc_func = "hclust", hc_method = "complete", method = "gap_stat")+
  labs(subtitle = "Gap statistic method")
# k=2
# Par la suite, nous prendrons arbitrairement k=4.

#Dendogramme pour le nombre de classe sélectionné.
par(mfrow=c(1,1))
fviz_dend(agnes.res, 4, cex=0.6, main = "Dendogramme de la partition en 4 classes", xlab = 'Individus')

# Graphique de classification associé
s.class(data,as.factor(cutree(agnes.res, k=4)), col = brewer.pal(5, "Set1"), sub = "Axes 1 et 2")
clusplot(data, cutree(agnes.res, k=4))



    # Optionnel : Fonction HCPC() pour CAH

    # Comparaison CAH sur données brutes et sur ACP
# CAH sur données centrées réduites
hcpc1.res <- HCPC(data.frame(scale(na.omit(eaux2018)[,3:11])), nb.clust=-1)
# CAH sur ACP
hcpc2.res <- HCPC(dudi.pca(data, scannf=F, nf=9)$li, nb.clust=-1)
cor(dist(hcpc1.res$data.clust),dist(data))
#0.997 : entre CAH sur données brutes et données brutes.
cor(dist(hcpc2.res$data.clust),dist(data))
#0.997 : entre CAH sur ACP et données brutes.
cor(dist(hcpc1.res$data.clust),dist(hcpc2.res$data.clust))
#1 : entre CAH sur données brutes et CAH sur ACP.

    # Comparaison CAH sur ACP consolidée et non-consolidée
acp = PCA(data, ncp = 5, graph = F)
res = HCPC(acp, consol = T) #une consolidation kmeans est réalisée
res$data.clust
res2 = HCPC(acp, consol = F) #pas de consolidation kmeans
res2$data.clust
cor(dist(res$data.clust),dist(data))
#0.997 : entre CAH sur données consolidées et données brutes.
cor(dist(res2$data.clust),dist(data))
#0.995 : entre CAH sur données non-consolidées et données brutes.
cor(dist(res$data.clust),dist(res2$data.clust))
# Corrélation de 0.995 entre la CAH sur ACP consolidée par kmeans et ACP non consolidée.



    # Classification descendante hiérarchique (DIANA - DIvisive ANAlysing)

diana.res <- diana(data)

#Sauts d'inertie
par(mfrow=c(1,1))
inertie <- sort(diana.res$height, decreasing = TRUE)
plot(inertie[1:20], type = "s", xlab = "Nombre de classes", ylab = "Inertie", main= "Graphique des sauts d'inertie")
points(3, inertie[3], col = c("green3","red3"), cex = 2, lwd = 3)
#k=3 se distingue des autres

#Calcul de l'inertie relative
best.cutree(diana.res, graph = TRUE, xlab = "Nombre de classes", ylab = "Inertie relative")
#La classification k=3 est devant toutes les autres.

# average silhouette ("silhouette"), within cluster sums of squares ("wss"), gap statistics ("gap_stat").
fviz_nbclust(data, hcut, method = "silhouette")+
  labs(subtitle = "Silhouette method") 
# k=3
fviz_nbclust(data, hcut, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method") 
# k=3
fviz_nbclust(data, hcut, method = "gap_stat")+
  labs(subtitle = "Gap statistic method")
# k=10
# Par la suite, nous prendrons arbitrairement k=3.

#Dendogramme pour le nombre de classe sélectionné.
par(mfrow=c(1,1))
fviz_dend(diana.res, 3, cex=0.6, main = "Dendogramme de la partition en 3 classes", xlab = 'Individus')

# Graphique de classification associé
s.class(data,as.factor(cutree(agnes.res, k=3)), col = brewer.pal(5, "Set1"), sub = "Axes 1 et 2")
clusplot(data, cutree(agnes.res, k=3))




    # Optionnel : Double CAH

heatmap(data, scale='none')

# Heatmap avec les variables en colonne et individus en ligne
Heatmap(data, name = "", column_title = "Variables", row_title = "Individus",
          row_names_gp = gpar(fontsize = 7))

# Le Biclustering ne nous permet pas dans notre cas d'arriver à une conclusion,
# quant à l'existence d'une structure sous-jacente entre les lignes et les colonnes du tableau de données.







    #################### Mesure de validité interne et externe ####################


      # Mesures interne (silhouette, Dunn et corrélation cophénétique)

data <- scale(na.omit(eaux2018)[,3:11])

    #kmeans
  #Silhouette (bon : proche de 1)
km.res1 <- eclust(data, "kmeans", k=2)
fviz_silhouette(km.res1, palette = "jco", ggtheme = theme_classic())
sil <- km.res1$silinfo$widths[, 1:3]
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]
# L'individu 75 est entre deux classes.
km.res1$silinfo$avg.width #Silhouette moyenne : 0.452
  #Indice de Dunn (indice à maximiser)
km_stats <- cluster.stats(dist(data), km.res1$cluster)
km_stats$dunn #L'indice de Dunn est de 0.278
  # Connectivité (indice à minimiser)
km_co <- connectivity(dist(data), km.res1$cluster)
km_co
#L'indice de connectivité de la classification est de 7.97

    #pam
  #Silhouette (bon : proche de 1)
pam.res1 <- eclust(data, "pam", k=3)
fviz_silhouette(pam.res1, palette = "jco", ggtheme = theme_classic())
sil <- pam.res1$silinfo$widths[, 1:3]
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]
# Les individus 75 et 12 sont entre les classes 1 et 2.
pam.res1$silinfo$avg.width #Silhouette moyenne : 0.516
  #Indice de Dunn (indice à maximiser)
pam_stats <- cluster.stats(dist(data), pam.res1$cluster)
pam_stats$dunn #L'indice de Dunn est de 0.320
  # Connectivité (indice à minimiser)
pam_co <- connectivity(dist(data), pam.res1$clustering)
pam_co
#L'indice de connectivité de la classification est de 11.71

    #clara
  #Silhouette (bon : proche de 1)
clara.res1 <- eclust(data, "clara", k=3)
fviz_silhouette(clara.res1, palette = "jco", ggtheme = theme_classic())
sil <- clara.res1$silinfo$widths[, 1:3]
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]
# Les individus 75 et 12 sont entre les classes 1 et 2.
clara.res1$silinfo$avg.width #Silhouette moyenne : 0.516
  #Indice de Dunn (indice à maximiser)
clara_stats <- cluster.stats(dist(data), clara.res1$cluster)
clara_stats$dunn #L'indice de Dunn est de 0.320
  # Connectivité (indice à minimiser)
clara_co <- connectivity(dist(data), clara.res1$clustering)
clara_co
#L'indice de connectivité de la classification est de 11.71

    #fanny
  #Silhouette (bon : proche de 1)
fanny.res1 <- eclust(data, "fanny", k=2)
fviz_silhouette(fanny.res1, palette = "jco", ggtheme = theme_classic())
sil <- fanny.res1$silinfo$widths[, 1:3]
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]
# Les individus 15, 19, 77, 20, 12, 58, 75, 88, 33 et 64 sont entre deux classes.
fanny.res1$silinfo$avg.width #Silhouette moyenne : 0.419
  #Indice de Dunn (indice à maximiser)
fanny_stats <- cluster.stats(dist(data), fanny.res1$cluster)
fanny_stats$dunn #L'indice de Dunn est de 0.154
  # Connectivité (indice à minimiser)
fanny_co <- connectivity(dist(data), fanny.res1$clustering)
fanny_co
#L'indice de connectivité de la classification est de 17.09

    #agnes
  #Silhouette (bon : proche de 1)
agnes.res1 <- eclust(data, FUNcluster="hclust", hc_method="complete", k=4)
fviz_silhouette(agnes.res1, palette = "jco", ggtheme = theme_classic())
sil <- agnes.res1$silinfo$widths[, 1:3]
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]
# Les individus 61 et 15 sont entre les classes 1 et 2.
# L'individu 74 entre les classes 1 et 4.
agnes.res1$silinfo$avg.width #Silhouette moyenne : 0.493
  #Indice de Dunn (indice à maximiser)
agnes_stats <- cluster.stats(dist(data), agnes.res1$cluster)
agnes_stats$dunn #L'indice de Dunn est de 0.418
  # Connectivité (indice à minimiser)
agnes_co <- connectivity(dist(data), agnes.res1$cluster)
agnes_co
#L'indice de connectivité de la classification est de 18.23
  # Corrélation cophénétique (bon : proche de 1, et au-dessus de 0.75)
coph_agnes <- cophenetic(agnes.res1)
cc_agnes <- cor(dist(data), coph_agnes)
cc_agnes
#La corrélation cophénétique est bonne avec une valeurs de 0.872

    #diana
  #Silhouette (bon : proche de 1)
diana.res1 <- eclust(data, "diana", k=3)
fviz_silhouette(diana.res1, palette = "jco", ggtheme = theme_classic())
sil <- diana.res1$silinfo$widths[, 1:3]
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]
# L'individu 77 est entre les classes 1 et 2.
diana.res1$silinfo$avg.width #Silhouette moyenne : 0.410
  #Indice de Dunn (indice à maximiser)
diana_stats <- cluster.stats(dist(data), diana.res1$cluster)
diana_stats$dunn #L'indice de Dunn est de 0.187
  # Connectivité (indice à minimiser)
diana_co <- connectivity(dist(data), diana.res1$cluster)
diana_co
#L'indice de connectivité de la classification est de 16.16
  # Corrélation cophénétique (bon : proche de 1, et au-dessus de 0.75)
coph_diana <- cophenetic(diana.res1)
cc_diana <- cor(dist(data), coph_diana)
cc_diana
#La corrélation cophénétique est bonne avec une valeurs de 0.840

t2 <- round(data.frame(cbind(rbind(km_stats$dunn,km_co,km.res1$silinfo$avg.width), rbind(pam_stats$dunn,pam_co,pam.res1$silinfo$avg.width),
                             rbind(clara_stats$dunn,clara_co,clara.res1$silinfo$avg.width), rbind(fanny_stats$dunn,fanny_co,fanny.res1$silinfo$avg.width),
                             rbind(agnes_stats$dunn,agnes_co,agnes.res1$silinfo$avg.width), rbind(diana_stats$dunn,diana_co,diana.res1$silinfo$avg.width))),digits=3)
rownames(t2) <- c('Incide de Dunn','Indice de connectivité','Silhouette')
colnames(t2) <- c('kmeans', 'pam', 'clara', 'fanny', 'agnes', 'diana')
stat_table <- formattable(t2, list(area(row=1) ~ color_tile('transparent','lightgreen'),
                                   area(row=2) ~ color_tile('lightblue', 'transparent'),
                                   area(row=3) ~ color_tile('transparent','lightgreen')))
stat_table

# Conclusion validité interne :
# - Selon l'indice de Dunn, agnes semble être meilleure que les autres.
# - Selon l'indice de connectivité, kmeans semble être meilleure que les autres.
# - Slon les silhouettes, pam et clara semble être meilleure que les autres.
# - Nous remarquons que Pam et clara ont de bons indices (ils n'ont pas de mauvais indices).




        # Mesures externes (Indice de Rand ajusté, Meila's VI et Kappa de Cohen)

data <- data.frame(cbind(scale(na.omit(eaux2018[,2:12])[2:10]),na.omit(eaux2018[,2:12])$Nature))
colnames(data)[10]<-'Nature'
head(data)
  
  #Entre Nature et kmeans
table(eaux[names(km.res$cluster),2], km.res$cluster) #Matrice de confusion
# Les eaux plates sont majoritairement dans le groupe 2 (43 en G2 contre 1 en G1).
# Les eaux gazeuses sont majoritairement dans le groupe 1 (12 en G1 contre 6 en G2).
nat <- as.numeric(data$Nature)
clus_km <- cluster.stats(dist(data), nat, km.res$cluster, compareonly=T)
clus_km
#rand ajusté est à 0.570 et Meila's VI est à 0.612
  # Coefficient de Kappa de Cohen
km_ck <- cohen.kappa(table(eaux[names(km.res$cluster),2], km.res$cluster))
km_ck
# Le coefficient de Kappa est estimé à 0.7, ce qui est un 'bon' degré d'accord.

  
  #Entre Nature et pam
table(eaux[names(pam.res$cluster),2], pam.res$cluster) #Matrice de confusion
# Les eaux plates sont majoritairement dans le groupe 1 (39 en G1 contre 1 en G2 et 4 en G3).
# Les eaux gazeuses sont majoritairement dans le groupe 1 (13 en G2 contre 5 en G1 et 0 en G3).
nat <- as.numeric(data$Nature)
clus_pam <- cluster.stats(dist(data), nat, pam.res$cluster, compareonly=T)
clus_pam
#rand ajusté est 0.492 et Meila's VI est 0.773
  # Coefficient de Kappa de Cohen
pam_ck <- cohen.kappa(table(eaux[names(pam.res$cluster),2], pam.res$cluster))
pam_ck
# Le degré d'accord en les groupes 1 et 2 est 'bon'.
# Le degré d'accord en les groupes 1 et 3 est 'très mauvais'.
# Le degré d'accord en les groupes 2 et 3 est 'très mauvais'.
# Enfin, la moyenne du coefficient de Kappa de Cohen est de 0.04 ce qui est 'mauvais'.

  #Entre Nature et clara
table(eaux[names(clara.res$cluster),2], clara.res$cluster) #Matrice de confusion
#la matrice de confusion est idem à pam
nat <- as.numeric(data$Nature)
clus_clara <- cluster.stats(dist(data), nat, clara.res$cluster, compareonly=T)
clus_clara
#les résultats sont idem à pam
#rand ajusté est 0.492 et Meila's VI est 0.773
  # Coefficient de Kappa de Cohen
clara_ck <- cohen.kappa(table(eaux[names(clara.res$cluster),2], clara.res$cluster))
clara_ck
# Le degré d'accord en les groupes 1 et 2 est 'bon'.
# Le degré d'accord en les groupes 1 et 3 est 'très mauvais'.
# Le degré d'accord en les groupes 2 et 3 est 'très mauvais'.
# Enfin, la moyenne du coefficient de Kappa de Cohen est de 0.04 ce qui est 'mauvais'.

  
  #Entre Nature et fanny
table(eaux[names(fanny.res$cluster),2], fanny.res$cluster) #Matrice de confusion
# Les eaux plates sont majoritairement dans le groupe 1 (38 en G1 contre 6 en G2).
# Les eaux gazeuses sont majoritairement dans le groupe 2 (16 en G2 contre 2 en G1).
nat <- as.numeric(data$Nature)
clus_fanny <- cluster.stats(dist(data), nat, fanny.res$cluster, compareonly=T)
clus_fanny
#rand ajusté est 0.538 et Meila's VI est 0.720
  # Coefficient de Kappa de Cohen
fanny_ck <- cohen.kappa(cbind(fanny.res$cluster,ifelse(eaux[names(fanny.res$cluster),2]=='gaz',2,1)))
fanny_ck
# Le coefficient de Kappa est estimé à 0.71, ce qui est un 'bon' degré d'accord.

  #Entre Nature et agnes
table(eaux[names(cutree(agnes.res, k=4)),2], cutree(agnes.res, k=4)) #Matrice de confusion
# Les eaux plates sont très majoritairement présentes dans le groupe 1 (39 en G1 contre 1 en G2 et 4 en G3).
# Les eaux gazeuses sont majoritairement dans le groupe 1 (11 en G1 contre 6 en G2 et 1 en G4).
# Néanmoins, plat est majoritaire en G1 et G2, tandise que gaz est majoritaire en G2 et G4.
nat <- as.numeric(data$Nature)
clus_agnes <- cluster.stats(dist(data), nat, cutree(agnes.res, k=4), compareonly=T)
clus_agnes
#rand ajusté est 0.200 et Meila's VI est 1.00
  # Coefficient de Kappa de Cohen
agnes_ck <- cohen.kappa(table(eaux[names(cutree(agnes.res, k=4)),2], cutree(agnes.res, k=4)))
agnes_ck
# Le degré d'accord en les groupes 1 et 2 est 'médiocre'.
# Le degré d'accord en les groupes 1 et 3 est 'mauvais'.
# Le degré d'accord en les groupes 2 et 3 est 'mauvais'.
# Enfin, la moyenne du coefficient de Kappa de Cohen est de 0.19 ce qui est 'mauvais'.

  #Entre Nature et diana
table(eaux[names(cutree(agnes.res, k=4)),2], cutree(diana.res, k=3)) #Matrice de confusion (agnes et diana utilise les mêmes individus)
# Les eaux plates sont très majoritairement présentes dans le groupe 1 (42 en G1 contre 1 en G2 et 1 en G3).
# Les eaux gazeuses sont majoritairement dans le groupe 2 (11 en G2 contre 6 en G1 et 1 en G3).
nat <- as.numeric(data$Nature)
clus_diana <- cluster.stats(dist(data), nat, cutree(diana.res, k=3), compareonly=T)
clus_diana
#rand ajusté est 0.523 et Meila's VI est 0.763
 # Coefficient de Kappa de Cohen
diana_ck <- cohen.kappa(table(eaux[names(cutree(agnes.res, k=3)),2], cutree(diana.res, k=3)))
diana_ck
# Le degré d'accord en les groupes 1 et 2 est 'bon'.
# Le degré d'accord en les groupes 1 et 3 est 'mauvais'.
# Le degré d'accord en les groupes 2 et 3 est 'très mauvais'.
# Enfin, la moyenne du coefficient de Kappa de Cohen est de 0.24 ce qui est 'médiocre'.


  #Création d'un tableau contenant les indices de Rand ajustés et Meila's VI
t1 <-round(data.frame(cbind(rbind(as.numeric(clus_km[1]),as.numeric(clus_km[2]),km_ck$weighted.kappa), 
                            rbind(as.numeric(clus_pam[1]),as.numeric(clus_pam[2]),pam_ck$av.wt),
                            rbind(as.numeric(clus_clara[1]),as.numeric(clus_clara[2]),clara_ck$av.wt),
                            rbind(as.numeric(clus_fanny[1]),as.numeric(clus_fanny[2]),fanny_ck$weighted.kappa),
                            rbind(as.numeric(clus_agnes[1]),as.numeric(clus_agnes[2]),agnes_ck$av.wt),
                            rbind(as.numeric(clus_diana[1]),as.numeric(clus_diana[2]),diana_ck$av.wt))), digits=3)
rownames(t1) <- c('Rand ajusté', 'Meila\'s VI', 'Kappa de Cohen')
colnames(t1) <- c('kmeans', 'pam', 'clara', 'fanny', 'agnes', 'diana') 
clus_table <- formattable(t1, list(
  area(row = 1) ~ color_tile('transparent','lightgreen'),
  area(row = 2) ~ color_tile('lightblue','transparent'),
  area(row = 3) ~ color_tile('transparent','lightgreen')))
clus_table

# Conclusion de validité externe : 
# - Selon rand ajusté, kmeans semble meilleure que les autres.
# - Selon Meila's VI, kmeans semble meilleure que les autres.
# - Selon Kappa de Cohen, fanny semble meilleure que les autres.
# - Nous remarquons que kmeans est clairement le meilleur algorithme pour la validation externe.



    # Choix de l'algorithme après validation 'manuelle'
x <- matrix(0, nrow=1, ncol=6)
x <- round((rank(data.frame(stat_table)[1,], ties.method='max') + (6-rank(data.frame(stat_table)[2,], ties.method='max')) +
        rank(data.frame(stat_table)[3,], ties.method='max') +
        rank(data.frame(clus_table)[1,], ties.method='max') + (6-rank(data.frame(clus_table)[2,], ties.method='max')) +
        rank(data.frame(clus_table)[3,], ties.method='max'))/6)
x <- formattable(data.frame(t(x)))
x
# -> D'après x qui nous donne le rang moyen (arrondi) obtenu pour les indices internes et externes présent dans clus_table et stat_table,
# kmeans arrive premier avec une constance dans la qualité de ses indices.



      # Choix de l'algorithme avec clValid

  # Mesures internes de clValid
intern <- clValid(data, 2:8, clMethods = c('kmeans','fanny','pam','clara', 'agnes', 'diana'), 
                 validation = c('internal'), method = 'complete')
summary(intern)

op <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(intern, legend = FALSE)
plot(nClusters(intern), measures(intern, "Dunn")[, , 1], type = "n", axes = F, 
     xlab = "", ylab = "")
legend("center", clusterMethods(intern), col = 1:6,lty = 1:3,pch = paste(1:6))
par(op)

# L'indice de connectivité nous propose kmeans avec k=2, c'est bien le choix que nous avons fait précédement.
# Nous remarquons que sur les trois indices, deux nous proposent kmeans, ce qui renforce notre choix effectué précédement.


  #Mesures de stabilité
stab <- clValid(data, 2:6, clMethods = c('kmeans','fanny','pam','clara', 'agnes', 'diana'), 
                 validation = 'stability', method = 'complete')
summary(stab)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(stab, measure = c("APN", "AD", "ADM"), legend = FALSE)
plot(nClusters(stab), measures(stab, "APN")[, , 1],type = "n", axes = F, 
     xlab = "", ylab = "")
legend("center", clusterMethods(stab), col = 1:6, lty = 1:3, pch = paste(1:6))
par(op)

# Les indices de stabilité nous orienterais plus vers clara avec k=3 (APN et ADM).


# CONCLUSION : notre choix s'étant fait par étape (choix de k, puis validation interne et externe),
# il paraît normal d'obtenir des résultats différents avec clValid qui lui calcul toutes les possiblités en même temps.




########## Graphique d'interprétation des résultats ##########

  # Starplot des 2 clusters obtenus

data <- data.frame(scale(na.omit(eaux2018)[,3:11]))

# Starplot élaboré des moyennes des deux clusters
star_moy <- data.frame(rbind(rep(3,9) , rep(-1,9) , round(colMeans(data[km.res$cluster==1,]),3), round(colMeans(data[km.res$cluster==2,]),3)))
star_moy
radarchart(star_moy, axistype=1, pcol=c('#987778','#166661'), pfcol=c('#987778','#166661'), plwd=2 ,
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(-1,3,1), cglwd=0.8,
            vlcex=0.8, title='Radar des moyennes pour chaque variable des 2 clusters')
legend(x='bottomleft', legend=c('Cluster 1   ','Cluster 2   '), fill=c('#987778','#166661'))

# Starplot élaboré des médianes des deux clusters
star_me <- matrix(0, nrow=4, ncol=9)
for (i in 1:9){
  star_me[1,i] <- 3
  star_me[2,i] <- -1
  star_me[3,i] <- round(median(data[km.res$cluster==1,i]),3)
  star_me[4,i] <- round(median(data[km.res$cluster==2,i]),3)
}
star_me <- data.frame(star_me)
colnames(star_me)<-colnames(data)
star_me
radarchart(star_me, axistype=1, pcol=c('#987778','#166661'), pfcol=c('#987778','#166661'), plwd=2 ,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(-1,3,1), cglwd=0.8,
           vlcex=0.8, title='Radar des médianes pour chaque variable des 2 clusters')
legend(x='bottomleft', legend=c('Cluster 1   ','Cluster 2   '), fill=c('#987778','#166661'))

# Starplot élaboré des écart-types des deux clusters
star_sd <- matrix(0, nrow=4, ncol=9)
for (i in 1:9){
  star_sd[1,i] <- 2
  star_sd[2,i] <- 0
  star_sd[3,i] <- round(sd(data[km.res$cluster==1,i]),3)
  star_sd[4,i] <- round(sd(data[km.res$cluster==2,i]),3)
}
star_sd <- data.frame(star_sd)
colnames(star_sd)<-colnames(data)
star_sd
radarchart(star_sd, axistype=1, pcol=c('#987778','#166661'), pfcol=c('#987778','#166661'), plwd=2 ,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,2,0.5), cglwd=0.8,
           vlcex=0.8, title='Radar des écart-types pour chaque variable des 2 clusters')
legend(x='bottomleft', legend=c('Cluster 1   ','Cluster 2   '), fill=c('#987778','#166661'))


  # Diagrammes en boîte des deux clusters
par(mfrow=c(1,2))
boxplot(data[km.res$cluster==1,], main=paste('Cluster 1 - ', nrow(data[km.res$cluster==1,]), 'indiv')) #cluster 1
boxplot(data[km.res$cluster==2,], main=paste('Cluster 2 - ', nrow(data[km.res$cluster==2,]), 'indiv')) #cluster 2
par(mfrow=c(1,1))




    # Matrice de confusion entre les 2 clusters kmeans / 3 groupes hypothétiques (plat fr/gaz fr/plat mar).

data <- data.frame(cbind(scale(na.omit(eaux2018[,2:12])[2:10]),na.omit(eaux2018[,2:12])[,c(1,11)]))
colnames(data)[10]<-'Nature'
colnames(data)[11]<-'Pays'
head(data)

table(km.res$cluster, ifelse(data$Pays==1, ifelse(data$Nature==0, 2, 1), 3)) #toutes les eaux marocaines sont plates
# Nous remarquons que les eaux plates qu'elles soient franaises ou marocaines, sont majoritairement dans le cluster 2,
# tandise que les eaux gazeuses sont majoritairement dans le cluster 1 mais la répartition est moins nette.








######################## Partie 2 - Classification supervisée ######################## 


# Cette seconde partie vise à mettre en place des techniques de classification supervisée pour prévoir la nature de l'eau.
# Pour ce faire, nous utiliserons les méthodes suivantes:
# - Méthode des k plus proche voisin
# - Analyse discriminante
# - Arbre de décision binaire


data <- data.frame(scale(na.omit(eaux2018)[,3:11]))

dim(data) #62 individus et 9 variables



    #Réalisons une fonction d'erreur de prédiction:
#Taux d'erreur de classement (nombre d'individus mal affectés sur le nombre total d'individus).
fct.err <- function(y, yhat){
  mc<-table(y,yhat)
  errare <- sum(1*(yhat!=y))
  #Calcul la valeur d'erreur
  return(errare/sum(mc))
}






        # Méthode des k plus proche voisin (kPPV) ou k nearest neighborhood (kNN)


    # Méthode d'échantillon/test

data <- data.frame(cbind(scale(na.omit(eaux2018[,2:12])[2:10]),na.omit(eaux2018[,2:12])$Nature))
colnames(data)[10]<-'Nature'
head(data)
#Séparons le dataset en train et test (80% vs 20% respectivement)
set.seed(13)
size_spl <- round(0.8*nrow(data))
train_ind <- sample.int(nrow(data), size=size_spl)
train <- data[train_ind,]
test <- data[-train_ind,]
all.err<-matrix(0,ncol=9,nrow=1)
for (j in 1:9){
    #apprendre le modèle à tous les individus
    knn_model<-knnreg(Nature ~ ., data=train, k=j)
    #appliquer le modèle sur le set de test
    pred<-round(predict(knn_model,newdata=test[,1:9]))
    #Calcul des erreurs par la fonction d'erreur
    all.err[1,j]<-fct.err(test$Nature, pred)
}

#Vecteur des erreurs pour chaque k plus proche voisin
print(round(all.err, 4))
et_KNN <- min(round(all.err, 4))
et_KNN
# Nous obtenons une erreur de 0.083 par méthode de train/test pour k={3,4,5,7} plus proches voisins.
# C'est un bon résultat quant à la prédiction du modèle kNN.


    # Validation croisée stratifiée

#Déterminer le numéro de bloc de chaque individu
n <- nrow(data) #nombre d'observations
k <- 7 #pour 8-validations croisées
taille <- n%/%k #déterminer la taille de chaque bloc
set.seed(13) #pour obtenir la même séquence
alea<-runif(n) #générer une colonne de valeurs aléatoires
rang<-rank(alea) #associer un individu à chaque rang
bloc<-(rang-1)%/%taille +1 #associer chaque individu à un numéro de bloc
bloc<-factor(bloc) #transformer en factor
summary(bloc) #vérification des classes créées

#Lancer la validation croisée stratifiée pour Nature
data <- data.frame(cbind(scale(na.omit(eaux2018[,2:12])[2:10]),na.omit(eaux2018[,2:12])$Nature))
colnames(data)[10]<-'Nature'
head(data)
all.err<-matrix(0,ncol=9,nrow=8)
for (j in 1:9){
for (i in 1:8){
  #apprendre le modèle à tous les individus sauf le blocs i
  knn_model<-knnreg(Nature ~ ., data=data[bloc!=i,], k=j)
  #appliquer le modèle sur le bloc i
  pred<-round(predict(knn_model,newdata=data[bloc==i,]))
  #Calcul des erreurs par la fonction d'erreur
  all.err[i,j]<-fct.err(data[bloc==i,10], pred)
}}

#Matrice des erreurs
print(round(all.err, 4))
#Moyenne des erreurs pour chaque k plus proche voisin
err.cv<-round(colMeans(all.err), 4)
err.cv
vcs_KNN <- min(err.cv)
vcs_KNN
#Nous obtenons une erreur moyenne de 0.094 par validation croisée stratifiée avec k={3,5} plus proches voisins.
# C'est un bon résultat quant à la prédiction du modèle kNN.



    # Leave or out (LOO)

head(data)
n <- nrow(data)
all.err<-matrix(0,nrow=n,ncol=9)
for (j in 1:9){
for (i in 1:n){
  # Apprendre le modèle à tous les individus étrangers au bloc i
  discr <- knnreg(Nature ~ ., data[-i,], k=j)
  # Calcul de l'erreur
  pred <- round(predict(discr,data[i,-10]))
  # Taux d'erreur de kNN
  all.err[i,j] <- fct.err(data[i,10], pred)
}}
all.err
colMeans(all.err)
loo_KNN <- min(colMeans(all.err))
loo_KNN # 0.097



    # Technique du bootstrap

data <- data.frame(cbind(scale(na.omit(eaux2018[,2:12])[2:10]),na.omit(eaux2018[,2:12])$Nature))
colnames(data)[10]<-'Nature'
head(data)

#Réalisons un bootstrap (rééchantillonage)
boot_KNN <- matrix(0, ncol=9)
for (j in 1:9){
  set.seed(13) # Pour obtenir le même résultat
  boot_error_rates <- sapply(seq_len(1000), function(b) {
    training <- sample(seq_along(data[,10]), replace = TRUE)
    train_out <- knnreg(data[training,-10], data[training,10], k=j)
    classifications <- predict(train_out, data[,-10])
    mean(classifications != data[,10])
  })
  boot_KNN[1,j] <- mean(boot_error_rates)
}
round(boot_KNN, 4)
bs_KNN <- min(round(boot_KNN, 4))
bs_KNN
# Nous obtenons une erreur de minimal de 0.059 pour k=1 et qui augmente plus k augmente.




    #Représentation graphique des prédiction par rapport aux valeurs réelles
plot(rownames(test), round(predict(knnreg(Nature ~ .,data=train, k=3),newdata=test[,1:9])), xlab="y", ylab=expression(hat(y)))
points(rownames(test), test$Nature, col='red', pch=5)
# L'individu 70 est le seul individu mal prévu (plat au lieu de gaz).






        # Analyse discriminante

  # Test pour l'optimalité de l'analyse discriminante
data <- na.omit(eaux2018)[,3:11]
head(data)
dat <- cbind(scale(data),na.omit(eaux2018)[,2])
colnames(dat)[10] <- 'Nature'
head(dat)
# Test de colinéarité des variables (indépendance entre les variables)
chisq.test(data[,2],data[,3])$p.value # On rejette H0 (hypothèse de non-colinéarité)
# Test d'homoscédasticité (homogénéité des variances)
bartlett.test(data) # On rejette H0 (hypothèse d'homoscédasticité)
# "qda"
# Test de normalité
for (i in 1:9){
  print(shapiro.test(data[,i])$p.value) # On rejette H0 (hypothèse de normalité)
}
mshapiro.test(t(data)) # On rejette H0 (hypothèse de multinormalité)




    # Méthode d'échantillon/test

data <- data.frame(cbind(scale(na.omit(eaux2018[,2:12])[2:10]),na.omit(eaux2018[,2:12])$Nature))
colnames(data)[10]<-'Nature'
data$Nature <- ifelse(data$Nature==1, 'gaz', 'plat')
head(data)
#Séparons le dataset en train et test (80% vs 20% respectivement)
set.seed(13)
size_spl <- round(0.8*nrow(data))
train_ind <- sample.int(nrow(data), size=size_spl)
train <- data[train_ind,]
test <- data[-train_ind,]
#apprendre le modèle à tous les individus
discr_model<-lda(Nature ~ ., data=train)
#appliquer le modèle sur le set de test
pred<-predict(discr_model,newdata=test[,1:9])$class
#Taux d'erreur de l'analyse discriminante
et_AD <- fct.err(test$Nature, pred)
et_AD
#Nous obtenons une erreur généralement inférieur à 0.083 par méthode de train/test.
# C'est un bon résultat quant à la prédiction du modèle d'analyse discriminante.



    # Validation croisée stratifiée

#Déterminer le numéro de bloc de chaque individu
n <- nrow(data) #nombre d'observations
k <- 7 #pour 8-validations croisées
taille <- n%/%k #déterminer la taille de chaque bloc
set.seed(13) #pour obtenir la même séquence
alea<-runif(n) #générer une colonne de valeurs aléatoires
rang<-rank(alea) #associer un individu à chaque rang
bloc<-(rang-1)%/%taille +1 #associer chaque individu à un numéro de bloc
bloc<-factor(bloc) #transformer en factor
summary(bloc) #vérification des classes créées

#Lancer la validation croisée stratifiée pour Nature
data <- data.frame(cbind(scale(na.omit(eaux2018[,2:12])[2:10]),na.omit(eaux2018[,2:12])$Nature))
colnames(data)[10]<-'Nature'
data$Nature <- ifelse(data$Nature==1, 'gaz', 'plat')
head(data)
all.err<-matrix(0,nrow=8,ncol=1)
for (i in 1:8){
  #apprendre le modèle à tous les individus étrangers au bloc i
  discr_model<-lda(Nature ~ ., data=data[bloc!=i,])
  #appliquer le modèle
  pred<-predict(discr_model,newdata=data[bloc==i,])$class
  #Taux d'erreur de l'analyse discriminante
  all.err[i] <- fct.err(data[bloc==i,10], pred)
}
all.err
vcs_AD <- mean(all.err)
vcs_AD
#Nous obtenons une erreur moyenne de 0.109 par validation croisée stratifiée.
# C'est un bon résultat quant à la prédiction du modèle d'analyse discriminante.



    # Leave or out (LOO)

head(data)
n <- nrow(data)
all.err<-matrix(0,nrow=n,ncol=1)
for (i in 1:n){
  #apprendre le modèle à tous les individus étrangers au bloc i
  discr <- lda(Nature ~ ., data[-i,])
  # Calcul de l'erreur
  pred <- predict(discr,data[i,-10])$class
  #Taux d'erreur en analyse discriminante
  all.err[i,1] <- fct.err(data[i,10], pred)
}
all.err
loo_AD <- mean(all.err)
loo_AD # 0.081




    # Bootstrap

data <- data.frame(cbind(scale(na.omit(eaux2018[,2:12])[2:10]),na.omit(eaux2018[,2:12])$Nature))
colnames(data)[10]<-'Nature'
head(data)

set.seed(13) # Pour obtenir le même résultat
boot_error_rates <- sapply(seq_len(1000), function(b) {
  training <- sample(seq_along(data[,10]), replace = TRUE)
  train_out <- lda(Nature ~ ., data[training,])
  classifications <- predict(train_out, data[,-10])$class
  mean(classifications != data[,10])
})
bs_AD <- mean(boot_error_rates)
bs_AD
# L'estimation d'erreur de bootstrap est de 0.080
# C'est un bon résultat quant à la prédiction du modèle d'analyse discriminante.







        # Arbres de décision binaires (arbres de classement)


    # Méthode d'échantillon/test

data <- data.frame(cbind(scale(na.omit(eaux2018[,2:12])[2:10]),na.omit(eaux2018[,2:12])$Nature))
colnames(data)[10]<-'Nature'
data$Nature <- ifelse(data$Nature==1, 'gaz', 'plat')
head(data)

#Séparons le dataset en train et test (80% vs 20% respectivement)
set.seed(13)
size_spl <- round(0.8*nrow(data))
train_ind <- sample.int(nrow(data), size=size_spl)
train <- data[train_ind,]
test <- data[-train_ind,]
bina <- rpart(Nature ~ ., train)
plot(bina,compress=T,uniform=T,branch=0.4,margin=0.1)
text(bina)
# Calcul de l'erreur
pred <- predict(bina,test[,1:9])
B <- ifelse(max.col(pred)==1,'gaz','plat')
et_AB <- fct.err(test$Nature, B)
et_AB # 0.083

# Erreur méthode échantillon/test : KNN (0.083) = AD (0.083) = AB (0.083)



    # Validation croisée stratifiée

#Déterminer le numéro de bloc de chaque individu
n <- nrow(data) #nombre d'observations
k <- 7 #pour 8-validations croisées
taille <- n%/%k #déterminer la taille de chaque bloc
set.seed(13) #pour obtenir la même séquence
alea<-runif(n) #générer une colonne de valeurs aléatoires
rang<-rank(alea) #associer un individu à chaque rang
bloc<-(rang-1)%/%taille +1 #associer chaque individu à un numéro de bloc
bloc<-factor(bloc) #transformer en factor
summary(bloc) #vérification des classes créées

#Lancer la validation croisée stratifiée pour Nature
data <- data.frame(cbind(scale(na.omit(eaux2018[,2:12])[2:10]),na.omit(eaux2018[,2:12])$Nature))
colnames(data)[10]<-'Nature'
data$Nature <- ifelse(data$Nature==1, 'gaz', 'plat')
head(data)
all.err<-matrix(0,nrow=8,ncol=1)
for (i in 1:8){
  #apprendre le modèle à tous les individus étrangers au bloc i
  bina <- rpart(Nature ~ ., data[bloc!=i,])
  # Calcul de l'erreur
  pred <- predict(bina,data[bloc==i,-10])
  B <- ifelse(max.col(pred)==1,'gaz','plat')
  #Taux d'erreur de l'analyse discriminante
  all.err[i] <- fct.err(data[bloc==i,10], B)
}
all.err
vcs_AB <- mean(all.err)
vcs_AB
# Nous obtenons une erreur moyenne de 0.1094 par validation croisée stratifiée.
# C'est un bon résultat quant à la prédiction du modèle d'arbre binaire.


# Erreur validation croisée : AD (0.109) = AB (0.109) > KNN (0.094)



    # Leave or out (LOO)

head(data)
n <- nrow(data)
all.err<-matrix(0,nrow=n,ncol=1)
for (i in 1:n){
  #apprendre le modèle à tous les individus étrangers au bloc i
  bina <- rpart(Nature ~ ., data[-i,])
  # Calcul de l'erreur
  pred <- round(predict(bina,data[i,-10]))
  #Taux d'erreur de l'arbre binaire
  all.err[i,1] <- fct.err(data[i,10], pred)
}
all.err
loo_AB <- mean(all.err)
loo_AB # 0.113


# Erreur par LOO : AB (0.113) > KNN (0.097) > AD (0.081)




    # Bootstrap

data <- data.frame(cbind(scale(na.omit(eaux2018[,2:12])[2:10]),na.omit(eaux2018[,2:12])$Nature))
colnames(data)[10]<-'Nature'
head(data)

set.seed(13) # Pour obtenir le même résultat
boot_error_rates <- sapply(seq_len(1000), function(b) {
  training <- sample(seq_along(data[,10]), replace = TRUE)
  train_out <- rpart(Nature ~ ., data[training,])
  classifications <- round(predict(train_out, data[,-10]))
  mean(classifications != data[,10])
})
bs_AB <- mean(boot_error_rates)
bs_AB
# L'estimation d'erreur par bootstrap est de 0.102
# C'est un bon résultat quant à la prédiction du modèle d'arbre binaire.


# Erreur par bootstrap : AB (0.102) = AD (0.080) > KNN (0.059)




    # Table des taux d'erreurs

table_MVT <- round(matrix(c(et_KNN,vcs_KNN,loo_KNN,bs_KNN,et_AD,vcs_AD,loo_AD,bs_AD,et_AB,vcs_AB,loo_AB,bs_AB), ncol=3), 4)
rownames(table_MVT) <- c('Méthode échantillon/test','Validation croisée stratifiée', 'Validation croisée par LOO','Technique du bootstrap')
colnames(table_MVT) <- c('kNN','Analyse discriminante','Arbre binaire')
table_MVT <- formattable(data.frame(table_MVT))
table_MVT

# Conclusion : kNN est le meilleur algorithme de classification supervisée dans notre cas.




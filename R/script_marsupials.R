## codigos para carregar dados geneticos, funcionais, ambientais e de ocorrencia das especies
## de roedores em ec?tonos para analisar os fatores predizendo a variacao local
## na forma

## 236 individuos capturados
## 176 individuos sequenciados (Thaptomys nigrita do GenBank)
## 146 rodents

source("R/packages.R")

#------------------------------- #
#         OCCURRENCE 
#------------------------------- #

occur <- read.csv (here ("data","ocorrencia.csv"),h=T,sep=";",row.names=1)

#------------------------------- #
#        GENETICS 
# ------------------------------ #

# loading trees
load(here("output", "fitGTRopt.RodMar.Rdata"))
#load(here("output", "bs.GTR_phylo.Rdata"))

## obtain genetic distance from optimized tree
opt_gen_dist <- as.matrix(as.dist(cophenetic (fitGTRopt.RodMar$tree),diag=T))

# represeting frequency of distances
pdf(here("output", "hist_small.pdf"),height=4,width=4)

  hist(opt_gen_dist,xlab="Genetic distance", main="",breaks=50)
  abline (v=0.03,lwd=1,lty=2)
  abline (v=0.10,lwd=1,lty=2)
  abline (v=0.15,lwd=1,lty=2)

dev.off()


#------------------------------- #
#            MORPHOLOGY          
#------------------------------- #

## atributos dos individuos
ind_trait <- read.csv (here ("data",'ind_trait.csv'),h=T,sep=";")

## selecionar as esp?cies de sigmodontine
small_mammals <- unique (ind_trait$species)

## pegar somente os atributos das sp de sigmodontine
ind_trait <- ind_trait[which (ind_trait$species %in% small_mammals),]

## ------ subsetting

## remove size effect
atributo_sem_tamanho <-ind_trait [,3:7]/ind_trait[,9]
rownames(atributo_sem_tamanho) <- ind_trait[,1]

## pegar atributos dos ind sequenciados
atributo_sem_tamanho <- atributo_sem_tamanho [which(rownames(atributo_sem_tamanho) %in% 
                                                      rownames(as.matrix(opt_gen_dist))),]

## distancia entre individuos com atributos padronizados
atributo_sem_tamanho <- decostand (atributo_sem_tamanho[,c("ccauda","pataunha" ,"orelha")],
                                   method="standardize")

## distance
dist_morfo <- as.matrix(vegdist(atributo_sem_tamanho, method="euclidean",diag=T))

## pegar a distancia genetica entre ind de sp de sigmodontine
#DIST_GEN_SIGMODONTINE <- dist_gen [which(rownames(dist_gen) %in% ind_trait$individuo),which(colnames(dist_gen) %in% ind_trait$individuo)]
#DIST_GEN_SIGMODONTINE <- dist_gen [which(rownames(dist_gen) %in% ind_trait$individuo)]

#labels_sp <- substring(rownames(as.matrix(opt_gen_dist)),4,9)
#labels_sp <- gsub ("1","",labels_sp)
#labels_sp <- gsub ("0","",labels_sp)
#labels_sp <- gsub ("2","",labels_sp)
#labels_sp <- gsub ("3","",labels_sp)
#labels_sp <- gsub ("4","",labels_sp)
#labels_sp <- gsub ("5","",labels_sp)
#labels_sp <- gsub ("6","",labels_sp)
#labels_sp <- gsub ("7","",labels_sp)
#labels_sp <- gsub ("8","",labels_sp)
#labels_sp <- gsub ("9","",labels_sp)
#
## grafico de variacao intra vs inter
#unique_labels_sp <- unique(labels_sp)

# variacao intra especifica
#dist_intra <- lapply (unique_labels_sp, function (sp)
#	as.dist(opt_gen_dist [grep (sp, rownames(opt_gen_dist)),
#		grep (sp, rownames(opt_gen_dist))],diag=F,upper=F))

#media_intra <- unlist(lapply(dist_intra,mean))

# variacao interespecifica
#dist_inter <- lapply (unique_labels_sp, function (sp)
#  opt_gen_dist [grep (sp, rownames(opt_gen_dist)),
#			-grep (sp, rownames(opt_gen_dist))])

#dist_inter <- lapply(dist_inter, function (i)
#	apply(i,2,as.numeric))

#media_inter <- unlist(lapply(dist_inter,mean,na.rm=T))
#
#df_intra_inter <- rbind(data.frame (dist_genetica = media_intra, var="Intra"),
#	data.frame (dist_genetica = media_inter, var="Inter"))
#
#df_intra_inter$var <- as.factor(df_intra_inter$var)
#
#with(df_intra_inter, plot(var,dist_genetica,ylab= "Genetic distance",
#		xlab="Level of variation",pch=19,cex=0.5,col="gray70"))
#
### transformar a matriz completa em distancias
## dist_gen_sigm <- data.matrix (as.dist(opt_gen_dist,diag=F,upper=F))

#------------------------------- #
#            ENVIRONMENT          
#------------------------------- #

amb_data <- read.csv (here ("data","amb_data.csv"), h=T,sep=";")

cor (amb_data[,c(7:17)])

## distancia ambiental entre campo e floresta
#dist_amb <- vegdist (decostand (amb_data[,c("samam","herb","arbus","arvor","serap","abdoss","drh")],
#                                "standardize"), 
#                    "euclidean")
## transformar em matriz
#DIST_AMB <- as.matrix(dist_amb, diag=T,upper=T)
## variacao dentro de cada habitat
#amb_habitats <- unique (amb_data$Floresta.interface)
#var_intra_hab <- lapply (amb_habitats, function (hab)
#		DIST_AMB [which(amb_data$Floresta.interface == hab),
#				which(amb_data$Floresta.interface == hab)])
#
### variacao entre habitats
#var_inter_hab <- lapply (amb_habitats, function (hab)
#		DIST_AMB [which(amb_data$Floresta.interface == hab),
#				-which(amb_data$Floresta.interface == hab)])
#
#df_dif_habitat <- rbind(data.frame(distancia = as.vector(var_intra_hab[[1]]),comp="CxC"),
#	data.frame(distancia = as.vector(var_intra_hab[[2]]),comp="FxF"),
#	data.frame(distancia = as.vector(var_inter_hab[[1]]),comp="CxF"))
#
#with(df_dif_habitat, plot(as.factor(comp),distancia,ylab= "Environmental difference",
#		xlab="Level of variation",pch=19,cex=0.5,col="gray70"))
#

### distancia ambiental no ponto de ocorr?ncia dos ind
### pegar as ocorrencias das sp sequenciadas
occur <- occur [,which(colnames(occur) %in% rownames(as.matrix(opt_gen_dist)))]

## ambiente da ocorr?ncia da sp
amb_ocorr <- lapply (seq(1,ncol(occur)), function (especie)
  amb_data [occur[,especie] >= 1,c("samam","herb","arbus","arvor","serap","abdoss","drh")])

amb_ocorr_ind <- do.call(rbind, 
                         lapply (seq(1,length (amb_ocorr)), function (sp) 
                           
                           colMeans (amb_ocorr[[sp]]))
                         )

colnames(amb_ocorr_ind) <- colnames(amb_ocorr[[1]])

## todas as correlacoes < 0.7
cor(amb_ocorr_ind)
#
### distancia ambiental das variaveis padronizadas
amb_ocorr_ind_stand <- decostand(amb_ocorr_ind,method="standardize")
dist_amb <-as.matrix (vegdist(amb_ocorr_ind_stand,method="euclidean",diag=T))

#------------------------------- #
#        GEOGRAPHY 
# ------------------------------ #

## ambiente da ocorr?ncia da sp
geo_ocorr <- lapply (seq(1,ncol(occur)), function (especie)
  amb_data [occur[,especie] >= 1,c("lats","lonw")])

geo_ocorr_ind <- do.call(rbind, 
                         lapply (seq(1,length (geo_ocorr)), function (sp) 
                           
                           colMeans (geo_ocorr[[sp]]))
                         
                         )

# todas as correlacoes < 0.7
cor(amb_ocorr_ind)

## distancia geografica das variaveis padronizadas
dist_geo <- as.matrix(vegdist(geo_ocorr_ind,method="euclidean",diag=T))

#------------------------------- #
#        VARIANCE PARTITIONING
# ------------------------------ #

particao_var_ind <- varpart (as.numeric(dist_morfo),   
                             as.numeric(opt_gen_dist),
                             as.numeric(dist_amb), 
                             as.numeric(dist_geo),
                             add="cailliez")
plot(particao_var_ind)

# barplot

total <- data.frame (Genetics = particao_var_ind$part$indfract$Adj.R.square[1],
            Environment = particao_var_ind$part$indfract$Adj.R.square[2],
            Distance = particao_var_ind$part$indfract$Adj.R.square[3],
            GeneticsEnvironment = particao_var_ind$part$indfract$Adj.R.square[4],
            DistanceEnvironment = particao_var_ind$part$indfract$Adj.R.square[5],
            DistanceGenetics = particao_var_ind$part$indfract$Adj.R.square[6],
            DistanceEnvironmentGenetics = particao_var_ind$part$indfract$Adj.R.square[7],
            Residuals=particao_var_ind$part$indfract$Adj.R.square[8]
            )

##############################
#### particao da variancia entre individuos com varias dist geneticas
##############################
dist_threshold <- 0.03
dist_threshold2 <- 0.10
dist_threshold3 <- 0.15

# 0.24
particao_var_ind <- varpart ((dist_morfo[which(opt_gen_dist > dist_threshold3)]),
                             (opt_gen_dist[which(opt_gen_dist > dist_threshold3)]),
                             (dist_amb[which(opt_gen_dist > dist_threshold3)]), 
                             (dist_geo[which(opt_gen_dist > dist_threshold3)]),
                             add="cailliez")
plot(particao_var_ind)

fifteen <- data.frame (Genetics = particao_var_ind$part$indfract$Adj.R.square[1],
                       Environment = particao_var_ind$part$indfract$Adj.R.square[2],
                       Distance = particao_var_ind$part$indfract$Adj.R.square[3],
                       GeneticsEnvironment = particao_var_ind$part$indfract$Adj.R.square[4],
                       DistanceEnvironment = particao_var_ind$part$indfract$Adj.R.square[5],
                       DistanceGenetics = particao_var_ind$part$indfract$Adj.R.square[6],
                       DistanceEnvironmentGenetics = particao_var_ind$part$indfract$Adj.R.square[7],
                       Residuals=particao_var_ind$part$indfract$Adj.R.square[8]
)

# 0.14 to 0.24
particao_var_ind <- varpart ((dist_morfo[which(opt_gen_dist > dist_threshold2 & opt_gen_dist <= dist_threshold3)]),
                             (opt_gen_dist[which(opt_gen_dist > dist_threshold2 & opt_gen_dist <= dist_threshold3)]),
                             (dist_amb[which(opt_gen_dist > dist_threshold2 & opt_gen_dist <= dist_threshold3)]), 
                             (dist_geo[which(opt_gen_dist > dist_threshold2 & opt_gen_dist <= dist_threshold3)]),
                             add="cailliez")
plot(particao_var_ind)

ten_fifteen <- data.frame (Genetics = particao_var_ind$part$indfract$Adj.R.square[1],
                           Environment = particao_var_ind$part$indfract$Adj.R.square[2],
                           Distance = particao_var_ind$part$indfract$Adj.R.square[3],
                           GeneticsEnvironment = particao_var_ind$part$indfract$Adj.R.square[4],
                           DistanceEnvironment = particao_var_ind$part$indfract$Adj.R.square[5],
                           DistanceGenetics = particao_var_ind$part$indfract$Adj.R.square[6],
                           DistanceEnvironmentGenetics = particao_var_ind$part$indfract$Adj.R.square[7],
                           Residuals=particao_var_ind$part$indfract$Adj.R.square[8]
)

# 0.04 to 0.14
particao_var_ind <- varpart ((dist_morfo[which(opt_gen_dist > dist_threshold & opt_gen_dist <= dist_threshold2)]),
                             (opt_gen_dist[which(opt_gen_dist > dist_threshold & opt_gen_dist <= dist_threshold2)]),
                             (dist_amb[which(opt_gen_dist > dist_threshold & opt_gen_dist <= dist_threshold2)]), 
                             (dist_geo[which(opt_gen_dist > dist_threshold & opt_gen_dist <= dist_threshold2)]),
                             add="cailliez")
plot(particao_var_ind)

five_ten <- data.frame (Genetics = particao_var_ind$part$indfract$Adj.R.square[1],
                        Environment = particao_var_ind$part$indfract$Adj.R.square[2],
                        Distance = particao_var_ind$part$indfract$Adj.R.square[3],
                        GeneticsEnvironment = particao_var_ind$part$indfract$Adj.R.square[4],
                        DistanceEnvironment = particao_var_ind$part$indfract$Adj.R.square[5],
                        DistanceGenetics = particao_var_ind$part$indfract$Adj.R.square[6],
                        DistanceEnvironmentGenetics = particao_var_ind$part$indfract$Adj.R.square[7],
                        Residuals=particao_var_ind$part$indfract$Adj.R.square[8]
)


# intra
particao_var_ind <- varpart ((dist_morfo[which(opt_gen_dist <= dist_threshold)]),
                             (opt_gen_dist[which(opt_gen_dist <= dist_threshold)]),
                             (dist_amb[which(opt_gen_dist <= dist_threshold)]), 
                             (dist_geo[which(opt_gen_dist <= dist_threshold)]),
                             add="cailliez")
plot(particao_var_ind)

five <- data.frame (Genetics = particao_var_ind$part$indfract$Adj.R.square[1],
                    Environment = particao_var_ind$part$indfract$Adj.R.square[2],
                    Distance = particao_var_ind$part$indfract$Adj.R.square[3],
                    GeneticsEnvironment = particao_var_ind$part$indfract$Adj.R.square[4],
                    DistanceEnvironment = particao_var_ind$part$indfract$Adj.R.square[5],
                    DistanceGenetics = particao_var_ind$part$indfract$Adj.R.square[6],
                    DistanceEnvironmentGenetics = particao_var_ind$part$indfract$Adj.R.square[7],
                    Residuals=particao_var_ind$part$indfract$Adj.R.square[8]
)

## rbind

res <- rbind (total,
       fifteen,
       ten_fifteen,
       five_ten,
       five)

# cplot of components
pdf(here("output", "partitioning_result_small.pdf"),width = 5, height=5)
par(mfrow=c(1,1),mar=c(5, 4, 4, 6) + 0.1)
plot (seq (1,5),res[,1],type="l",
      pch=19,ylim=c(0,0.15),col="blue",
      axes=FALSE,xlab="", ylab="",
      lwd=2
      )

axis(2, ylim=c(0,1),col="black",las=1)
mtext("Adjusted R-square",side=2,line=3)
box()

# Allow a second plot on the same graph
par(new=TRUE)

plot (seq (1,5),res[,8],type="l",
      pch=19,ylim=c(0,1),col="gray50",
      axes=FALSE,xlab="", ylab="",
      lwd=2
)
## a little farther out (line=4) to make room for labels
mtext("Residuals",side=4,col="gray50",line=3) 
axis(4, ylim=c(0,1), col="gray50",col.axis="gray50",las=1)

lines(seq (1,5),res[,2],type="l",pch=19,col="green",
      lwd=2)
lines(seq (1,5),res[,3],type="l",pch=19,col="red",
      lwd=2)
lines(seq (1,5),res[,4],type="l",pch=19,col="yellow",
      lwd=2)
lines(seq (1,5),res[,5],type="l",pch=19,col="purple",
      lwd=2)
lines(seq (1,5),res[,6],type="l",pch=19,col="orange",
      lwd=2)
lines(seq (1,5),res[,7],type="l",pch=19,col="beige",
      lwd=2)

legend ("topright",c("Genetics (G)",
                     "Environment (E)",
                     "Distance (D)",
                     "G+E",
                     "D+E",
                     "D+G",
                     "D+G+E"),bty="n",
        lty=1,lwd=2, col=c("blue", 
                     "green",
                     "red",
                     "yellow",
                     "purple",
                     "orange",
                     "beige"
                     ),
        cex=0.75)

axis(1, xlim=c(0,5), col="black",col.axis="black",las=1)
mtext("Interval of genetic distance",side=1,col="black",line=3) 

dev.off()

## biplot of distance 
## most important components
## relacao dist morfologica e dist ambiental
pdf(here("output", "lm_morfo_genetic_amb_small.pdf"),heigh=4,width=8)

par(mfrow=c(1,2))
plot((dist_morfo) ~ (opt_gen_dist),pch=19,cex=0.5, 
     xlab="Genetic distance",ylab="Morphological disparity", col="gray30",
     main="Total genetic variation")
abline(lm (as.numeric(dist_morfo) ~ as.numeric(opt_gen_dist)),lwd=2,col="black")

summary(lm(dist_morfo ~ opt_gen_dist))

plot(dist_morfo[which(opt_gen_dist <= dist_threshold)] ~ 
       dist_amb[which(opt_gen_dist <= dist_threshold)],pch=19,cex=0.5, 
     xlab="Environmental distance",ylab="Morphological disparity", col="gray30",
     main="Genetic distance <=0.03")
abline((lm((dist_morfo[which(opt_gen_dist <= dist_threshold)]) ~ 
             (dist_amb[which(opt_gen_dist <= dist_threshold)]))),lwd=2,col="black")

summary (lm((dist_morfo[which(opt_gen_dist <= dist_threshold)]) ~ 
               (dist_amb[which(opt_gen_dist <= dist_threshold)])))

dev.off()


### alternative whay to analyze

intervals <- seq (0.03,range(opt_gen_dist)[2],  0.02) 

# intra

particao_var_ind <- lapply (intervals, function (intervals) 
  
                    varpart ((dist_morfo[opt_gen_dist <= intervals]),
                             (opt_gen_dist[opt_gen_dist <= intervals]),
                             (dist_amb[opt_gen_dist <= intervals]), 
                             (dist_geo[opt_gen_dist <= intervals]),
                             add="cailliez"))

# extract components

res_df <- lapply (particao_var_ind, function (i) 
  
        data.frame (Genetics = i$part$indfract$Adj.R.square[1],
                    Environment = i$part$indfract$Adj.R.square[2],
                    Distance = i$part$indfract$Adj.R.square[3],
                    GeneticsEnvironment = i$part$indfract$Adj.R.square[4],
                    DistanceEnvironment = i$part$indfract$Adj.R.square[5],
                    DistanceGenetics = i$part$indfract$Adj.R.square[6],
                    DistanceEnvironmentGenetics = i$part$indfract$Adj.R.square[7],
                    Residuals= i$part$indfract$Adj.R.square[8]
))


res_df <- do.call(rbind, res_df)

# plot

pdf(here("output", "continuous_genetics_small.pdf"),width = 5, height=5)

par(mfrow=c(1,1),mar=c(5, 4, 4, 6) + 0.1)
plot (intervals,res_df[,1],type="l",
      pch=19,ylim=c(0,0.15),col="blue",
      axes=FALSE,xlab="", ylab="",
      lwd=2
)

axis(2, ylim=c(0,1),col="black",las=1,cex.axis=0.85)
mtext("Adjusted R-square",side=2,line=3)
box()

# Allow a second plot on the same graph
par(new=TRUE)

plot (intervals,res_df[,8],type="l",
      pch=19,ylim=c(0,1),col="gray50",
      axes=FALSE,xlab="", ylab="",
      lwd=2
)
## a little farther out (line=4) to make room for labels
mtext("Residuals",side=4,col="gray50",line=3) 
axis(4, ylim=c(0,1), col="gray50",col.axis="gray50",las=1)

lines(intervals,res_df[,2],type="l",pch=19,col="green",
      lwd=2)
lines(intervals,res_df[,3],type="l",pch=19,col="red",
      lwd=2)
lines(intervals,res_df[,4],type="l",pch=19,col="yellow",
      lwd=2)
lines(intervals,res_df[,5],type="l",pch=19,col="purple",
      lwd=2)
lines(intervals,res_df[,6],type="l",pch=19,col="orange",
      lwd=2)
lines(intervals,res_df[,7],type="l",pch=19,col="beige",
      lwd=2)

legend ("topright",c("Genetics (G)",
                     "Environment (E)",
                     "Distance (D)",
                     "G+E",
                     "D+E",
                     "D+G",
                     "D+G+E"),bty="n",
        lty=1,lwd=2, col=c("blue", 
                           "green",
                           "red",
                           "yellow",
                           "purple",
                           "orange",
                           "beige"
        ),
        cex=0.75)

axis(1, xlim=c(0,0.7), cex.axis=0.85,col="black",col.axis="black",las=1)
mtext("Phylogenetic scale (Increasing genetic distance)",side=1,col="black",line=3) 

dev.off()


###############################################################################




################ DIFERENCA MORFOLOGICA E GENETICA ENTRE 
## INDIVIDUOS OCORRERNDO NO CAMPO, FLORESTA OU EM AMBOS

### ha diferenca da genetica dos ind de campos e florestas
amb_ocorr <- lapply (seq(1,ncol(occur)), function (especie)
  amb_data [occur[,especie] >= 1,"Floresta.interface"])

amb_ocorr_ind <- lapply(amb_ocorr, function (i) 
	paste(as.character(i)[1],as.character(i)[2],sep=""))
amb_ocorr_ind <- gsub("NA","",unlist(amb_ocorr_ind))
amb_ocorr_ind [which(amb_ocorr_ind == "ii")] <- "i"
amb_ocorr_ind [which(amb_ocorr_ind == "ff")] <- "f"

names(amb_ocorr_ind) <- colnames(occur)

matrix_factors <- as.dist(outer(as.numeric(as.factor(amb_ocorr_ind)), 
	rep(0,145), "+"),upper=F,diag=F)

## ent?o se somar
## 11 eh floresta
## 12 eh campo
## 13 eh ambos
## 14 

## individuos com dist genetica baixa

fator <- as.factor(matrix_factors [dist_gen_sigm <= 0.04])
levels(fator)[which(levels(fator) == "1")] <- "Floresta"
levels(fator)[which(levels(fator) == "2")] <- "Campo"
levels(fator)[which(levels(fator) == "3")] <- "Ambos"

## distancias entre todos os pares de ind
plot(as.factor(fator),
	dist_gen_sigm [dist_gen_sigm <= 0.04],
     pch=19,cex=0.5,col="gray30",xlab="Habitat",
	ylab="Dist?ncia gen?tica entre indiv?duos")

anova (lm(dist_gen_sigm [dist_gen_sigm <= 0.04] ~ fator))
summary(lm(dist_gen_sigm [dist_gen_sigm <= 0.04] ~ fator))
plot(lm(sqrt(dist_gen_sigm [dist_gen_sigm <= 0.04]) ~ fator))


plot(as.factor(fator),
	dist_morfo [dist_gen_sigm <= 0.04],
     pch=19,cex=0.5,col="gray30",xlab="Habitat",
	ylab="Dist?ncia morfol?gica entre indiv?duos")
anova (lm(dist_morfo [dist_gen_sigm <= 0.04] ~ fator))
summary(lm(dist_morfo [dist_gen_sigm <= 0.04] ~ fator))


################################ 
## ANALISES PARA COMUNIDADES  ##
################################

occur <- decostand (occur,"pa")
occur <- occur [rowSums(occur) >=2,]

## lista de comunidades com presencas
list_occur <- apply(occur,1,list)
list_occur <-lapply (seq(1,length (list_occur)), function (comm)
  list_occur[[comm]][[1]][which(list_occur[[comm]][[1]] ==1)])

## atributos dos ind presentes em cada transeccao
comm_morfo <- lapply(list_occur, function (i) 
  atributo_sem_tamanho [which(rownames(atributo_sem_tamanho) %in% names(i)),])

## remover sem dados
comm_morfo <- comm_morfo [which(unlist(lapply(comm_morfo,nrow))>1)]

## distancia morfo entre ind coexistentes
comm_morfo_dist <- lapply (comm_morfo, function (comm)
  vegdist (comm, "euclidean"))

## distancia media entre ind coexistentes
media_dist_comm <- lapply (comm_morfo_dist, function (comm)
  mean (comm))

## distancia maxima entre ind coexistentes
max_dist_comm <- lapply (comm_morfo_dist, function (comm)
  max (comm))

## distancia minima entre ind coexistentes
min_dist_comm <- lapply (comm_morfo_dist, function (comm)
  min (comm))

## distancia genetica entre ind coexistentes
comm_gen_dist <- lapply(list_occur, function (i) 
  dist_gen_sigm [which(attr (dist_gen_sigm,"Labels") %in% names(i))])

## distancia gen media entre ind coexistentes
media_gen_com <- lapply (comm_gen_dist, function (comm)
  mean (comm))

## distancia gen maxima entre ind coexistentes
max_gen_comm <- lapply (comm_gen_dist, function (comm)
  max (comm))

## distancia gen minima entre ind coexistentes
min_gen_comm <- lapply (comm_gen_dist, function (comm)
  min (comm))

## ambiente para o subconj de comun com +2 ind coexistentes
amb_subset <- amb_data [which(rowSums(occur) >=2),c(10,14,15:17)]
amb_subset <- decostand (amb_subset,"standardize")

## relacao
## dist genetica
par(mfrow=c(2,3))
lapply (seq (1,ncol (amb_subset)), function (i)
  plot(amb_subset[,i], unlist(media_gen_com),pch=19,main=colnames(amb_subset)[i]))
## dist morfologica
par(mfrow=c(2,3))
lapply (seq (1,ncol (amb_subset)), function (i)
  plot(amb_subset[,i], unlist(media_gen_com),pch=19,main=colnames(amb_subset)[i]))

## disparidade morfologica
require(geomorph)

## geomorph dataframe
geom_morfo <- lapply (comm_morfo, function (comm)
  geomorph.data.frame (macaco=as.matrix(comm)))

## disparidade morfologica
## carregar a funcao modificada
load("morph_disp.RData")
morfo_disp <- lapply(geom_morfo, function (comm)
  morph_disp (macaco~1, data=comm, print.progress = T,seed=4,groups = NULL,partial=F))

##
par(mfrow=c(2,3))
lapply (seq (1,ncol (amb_subset)), function (i)
  plot(amb_subset[,i][which(unlist(lapply(comm_morfo,nrow))>1)], 
	unlist(morfo_disp),pch=19,main=colnames(amb_subset)[i]))

################### HAPLOVECTORS

## dados geneticos
require(ape)

gen_data <- read.dna("fas_gen_data.fas", format = "fasta")
ind_sequenciados <- gsub(" ","", attr (gen_data,"dimnames")[[1]])

## atributos dos individuos
ind_trait <- read.csv ('ind_trait.csv',h=T,sep=";")

## selecionar as esp?cies de sigmodontine
sigmodontine <- unique (ind_trait$species)[-c(1,4,9,10,16)]

## pegar somente os atributos das sp de sigmodontine
ind_trait <- ind_trait[which (ind_trait$species %in% sigmodontine),]

## remover a linha 31, ind de herval, que foi o unico ind da UAL 1
ind_trait <- ind_trait[-which(ind_trait$individuo == "hl1oxna"),]

## seleciona os sigmodontine dentro dos dados geneticos
gen_data  <- gen_data [which (ind_sequenciados %in% as.character(ind_trait$individuo) == T),]

## ocorrencia
occur <- read.csv ("ocorrencia.csv",h=T,sep=";",row.names=1)

## tirar o ind de herval
occur <- occur [,-which(colnames(occur) == "hl1oxna")]

## pegar as ocorrencias das sp sequenciadas
occur <- occur [,which(colnames(occur) %in% dimnames(gen_data)[1][[1]])]

### 'derreter' a matriz de occ em um df
require(reshape)
melt_occur <- cbind(occur,id=rownames(occur))
melt_occur <- melt(melt_occur,id="id")

## criando um fator "campo - floresta" e um fator sitio

fatores <- data.frame (habitat=rep(c(rep("campo", 9*8),rep("floresta", 9*8)),ncol(occur)),
		UAL = rep(substring(rownames(occur),2,2),ncol(occur)),
            sitio =rep(rep (c("APA","ENC", "HER", "PRO", "SAN","BOA", "ASS", "LIV","TAI"),16),ncol(occur)))

## colar os fatores ao df de ocorrencia
melt_occur <- cbind(melt_occur, fatores)
## colar uma interacao entre fatores
melt_occur <- cbind (melt_occur, interacao=paste(melt_occur$sitio,melt_occur$habitat,
	melt_occur$UAL,sep=""))

## uma matriz de ocorrencia de ind por cada combinacao da interacao
## remover a primeira coluna que refere-se ao fator
presaus <- cast(interacao ~variable, data= melt_occur,fun.aggregate=max)
rownames(presaus) <- presaus$interacao
presaus <- presaus[,-1]
## transformando em PA
presaus [presaus >= 1] <- 1

## pegando sitios com ao menos 2 ind
presaus <- presaus [which(rowSums(presaus)>1),]
## corrigindo a ordem das colunas
presaus <- presaus [,c(1:98,100,99,101:145)]

require(GenVectors)
require(pegas)
x<-haplotype(gen_data)
## gerar haplovetores
# teste <- HaploDist(gen_data)
# teste1 <- haplotype(gen_data) 

lista_obj <- list(gen =gen_data, occ = data.matrix(presaus),envir = as.factor(c(rep("c",10), rep("f",9))))

## composicao de haplotipos
comp_haplo <- HaploVectors(as.list(lista_obj$gen), lista_obj$occ)
plot(comp_haplo$vectors [,1],comp_haplo$vectors [,2],
	col=as.numeric(as.factor(substring (names(comp_haplo$vectors [,2]),4,7))),
	xlim=c (-0.2,0.2),ylim=c(-0.15,0.15))

for (i in 1:length(names(comp_haplo$vectors[,1]))) {
	text (comp_haplo$vectors [i,1],comp_haplo$vectors [i,2],
		names(comp_haplo$vectors[,1])[i],cex=0.65)
	}

plot(comp_haplo$vectors [,1],comp_haplo$vectors [,18],
	col=as.numeric(as.factor(substring (names(comp_haplo$vectors [,3]),4,7))),
	xlim=c (-0.2,0.2),ylim=c(-0.15,0.15))

for (i in 1:length(names(comp_haplo$vectors[,1]))) {
	text (comp_haplo$vectors [i,1],comp_haplo$vectors [i,18],
		names(comp_haplo$vectors[,1])[i],cex=0.65)
	}



## ocorrencias e calculo de disparidade morfologica
list_occur <- apply(presaus,1,list)
list_occur <- lapply (list_occur, function (i)
	{names(i [[1]])<- colnames(presaus);i})

list_occur <-lapply (seq(1,length (list_occur)), function (comm)
  list_occur[[comm]][[1]][which(list_occur[[comm]][[1]] ==1)])

comm_morfo <- lapply(list_occur , function (i) 
  atributo_sem_tamanho [which(rownames(atributo_sem_tamanho) %in% names(i)),])

require(betapart)

require(SYNCSA)

comp_functional <- matrix.t(data.matrix(presaus), atributo_sem_tamanho)$matrix.T
dist_functional <- as.matrix (vegdist (comp_functional,"euclidean"),upper=T,diag=T)

###########
## disparidade morfologica
#require(geomorph)

## geomorph dataframe
#geom_morfo <- lapply (comm_morfo, function (comm)
#  geomorph.data.frame (macaco=as.matrix(comm)))
## disparidade morfologica
## carregar a funcao modificada
#load("morph_disp.RData")
#morfo_disp <- lapply(geom_morfo, function (comm)
#  morph_disp (macaco~1, data=comm, print.progress = T,seed=4,groups = NULL,partial=F))

### correlacao da disparidade morfologica com composicao haplotipica
corr_morfo_gen_flor <- lapply (seq (1,ncol(comp_haplo$vectors)), function (i) 
	cor (apply(dist_functional,1,mean)[which(substring (names(comp_haplo$vectors [,1]),4,7)=="flor")],
	comp_haplo$vectors[,i][which(substring (names(comp_haplo$vectors [,1]),4,7)=="flor")]))

corr_morfo_gen_campo <- lapply (seq (1,ncol(comp_haplo$vectors)), function (i) 
	cor (apply(dist_functional,1,mean)[which(substring (names(comp_haplo$vectors [,1]),4,7)=="camp")],
	comp_haplo$vectors[,i][which(substring (names(comp_haplo$vectors [,1]),4,7)=="camp")]))

plot(unlist(corr_morfo_gen_flor),
	pch=19,type="b",cex=1.4,col="gray60",xlab= "HaploVector",
	ylab= "Correla??o linear",ylim=c(-1,1),xaxt="n")
points(unlist(corr_morfo_gen_campo),
	pch=19,type="b",cex=1.2,col="red1")
axis(seq(1,18),at=seq(1,18,1),cex.axis=1,las=2)
legend ("topleft", legend = c("Floresta", "Campo"),pch=19,
	col=c("gray60","red1"),bty="n")

campo_floresta_fator <- as.factor (substring (names(comp_haplo$vectors [,1]),4,7))
levels (campo_floresta_fator)[which(levels (campo_floresta_fator) == "camp")] <- "Campo"
levels (campo_floresta_fator)[which(levels (campo_floresta_fator) == "flor")] <- "Floresta"

plot(campo_floresta_fator, comp_haplo$vectors[,5])
plot(campo_floresta_fator, apply(dist_functional,1,max))

##


# barplot
#barplot (c(particao_var_ind$part$indfract$Adj.R.square[1],
#           particao_var_ind$part$indfract$Adj.R.square[2],
#           particao_var_ind$part$indfract$Adj.R.square[3],
#           particao_var_ind$part$indfract$Adj.R.square[4],
#           particao_var_ind$part$indfract$Adj.R.square[5],
#           particao_var_ind$part$indfract$Adj.R.square[6],
#           particao_var_ind$part$indfract$Adj.R.square[7]
#), names.arg = c("Genetics (G)",
#                 "Environment (E)",
#                 "Distance (D)",
#                 "G+E",
#                 "D+E",
#                 "D+G",
#                 "D+G+E"),
#col = c("blue",
#        "green",
#        "red",
#        "yellow",
#        "purple",
#        "orange",
#        "gray"),las=2, ylab="Ajusted R-square",
#ylim=c(0,0.4),
#main="Genetic distance <=0.06")
#text(x = 6,0.35,
#     paste ("Residuals=",
#            round (particao_var_ind$part$indfract$Adj.R.square[8],3)))





## test fraction
X1 <- rda(dist_morfo[dist_gen_sigm <= dist_threshold] ~
            dist_amb[dist_gen_sigm <= dist_threshold] + 
            Condition (dist_gen_sigm[dist_gen_sigm <= dist_threshold]) + 
            Condition (dist_geo[dist_gen_sigm <= dist_threshold]))
anova(X1)


dev.off()

## test of portions

## test fraction
X1 <- rda(dist_morfo ~
            dist_gen_sigm + 
            (dist_amb) + 
            (dist_geo))
anova(X1)

X2 <- rda(dist_morfo[dist_gen_sigm <= dist_threshold] ~
            Condition (dist_gen_sigm[dist_gen_sigm <= dist_threshold]) + 
            dist_amb[dist_gen_sigm <= dist_threshold] + 
            Condition (dist_geo[dist_gen_sigm <= dist_threshold]))
anova(X2)

X3 <- rda(dist_morfo[dist_gen_sigm <= dist_threshold] ~
            Condition (dist_gen_sigm[dist_gen_sigm <= dist_threshold]) + 
            Condition (dist_amb[dist_gen_sigm <= dist_threshold]) + 
            dist_geo[dist_gen_sigm <= dist_threshold])
anova(X3)

X4 <- rda(dist_morfo[dist_gen_sigm <= dist_threshold] ~
            (dist_gen_sigm[dist_gen_sigm <= dist_threshold]) + 
            (dist_amb[dist_gen_sigm <= dist_threshold]) + 
            Condition(dist_geo[dist_gen_sigm <= dist_threshold]))
anova(X4)

X5 <- rda(dist_morfo[dist_gen_sigm <= dist_threshold] ~
            Condition(dist_gen_sigm[dist_gen_sigm <= dist_threshold]) + 
            (dist_amb[dist_gen_sigm <= dist_threshold]) + 
            (dist_geo[dist_gen_sigm <= dist_threshold]))
anova(X5)

X6 <- rda(dist_morfo[dist_gen_sigm <= dist_threshold] ~
            (dist_gen_sigm[dist_gen_sigm <= dist_threshold]) + 
            Condition(dist_amb[dist_gen_sigm <= dist_threshold]) + 
            (dist_geo[dist_gen_sigm <= dist_threshold]))
anova(X6)

X7 <- rda(dist_morfo[dist_gen_sigm <= dist_threshold] ~
            (dist_gen_sigm[dist_gen_sigm <= dist_threshold]) + 
            (dist_amb[dist_gen_sigm <= dist_threshold]) + 
            (dist_geo[dist_gen_sigm <= dist_threshold]))
anova(X7)



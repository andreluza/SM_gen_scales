# ------------------------- #
# Genetic analysis

source("R/packages.R")

## help here : https://www.molecularecologist.com/2016/02/26/quick-and-dirty-tree-building-in-r/
# and here https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html

gen_data <- read.dna(here ("data","fas_gen_data.fas"), format = "fasta")
ind_sequenciados <- gsub(" ","", attr (gen_data,"dimnames")[[1]])

#------------------------------- #
#            MORPHOLOGY          
#------------------------------- #
## seleciona com base na presenca nos dados de traits
## atributos dos individuos
ind_trait <- read.csv (here ("data",'ind_trait.csv'),h=T,sep=";")

# subset
gen_data  <- gen_data [which (ind_sequenciados %in% as.character(ind_trait$individuo) == T),]

# DNA INTO PHYDAT
gen_small <- phyDat(gen_data, type = "DNA", levels = NULL)
gen_small <- subset (gen_small, attr (gen_data,"dimnames")[[1]] %in% ind_trait$individuo)
#class(gen_rodents)

## NJ tree
gen_small_dnabin <- as.DNAbin (gen_small)
# genetic distance
dist_gen <- dist.dna(gen_small_dnabin, model="K80") # distance derived by Kimura (1980),
small_NJ  <- nj(dist_gen)

## bootstrapping the tree
myBoots <- boot.phylo(small_NJ, gen_small_dnabin, function(e) 
  root(nj(dist.dna(e, model = "K80")),1))

plot(small_NJ, show.tip=T, edge.width=2)
title("NJ tree + bootstrap values")
tiplabels(frame="none", pch=20, col=rgb(0,0,0.1,alpha=0.1), cex=3, fg="transparent")
axisPhylo()
nodelabels(myBoots, cex=.6)

## select one evol model
#mt_rod_mars <- modelTest(gen_small)
#save(mt_rod_mars,file = here ("output", "test_model_rod_mars.RData"))

load(here ("output", "test_model_rod_mars.RData"))
env <- attr(mt_rod_mars, "env")
ls(envir=env)

round(mt_rod_mars$AICw,4)

## test diff models

(fitGTR <- eval(get("GTR+G+I", env), env))

## get the best and optimize tree
fitGTRopt.RodMar <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                       rearrangement = "stochastic", control = pml.control(trace = 0))
save(fitGTRopt.RodMar,file=here("output","fitGTRopt.RodMar.RData"))

# bootstrap to have branch support
nsamples<- 1000
bs.GTR.RodMar <- bootstrap.pml(fitGTRopt.RodMar, bs= nsamples, optNni=TRUE,
                        control = pml.control(trace = 0),
                        optGamma=TRUE, optInv=TRUE, model="GTR")

save(bs.GTR.RodMar,file=here("output","bs.GTR.RodMar.RData"))

## cons tree
#cons <- consensus(bs.GTR, p = 0.7, check.labels = TRUE)

load(here("output","bs.GTR.RodMar.RData"))
load(here("output","fitGTRopt.RodMar.RData"))
# plot opt tree and branch support

pdf(here("output","phyGTRFanRodMar_noTips.pdf"))
par(mfrow=c(1,1))
plotBS (midpoint(fitGTRopt.RodMar$tree),
        bs.GTR.RodMar,
        p = 5, 
        type="fan",
        cex=0.5,
        show.tip.label = F,
        edge.width = 1)
axisPhylo(side=2)

dev.off()

pdf(here("output","phyGTRRodMar.pdf"))
plotBS (midpoint(fitGTRopt.RodMar$tree),
        bs.GTR.RodMar,
        p = 5, 
        type="p",
        cex=0.5,
        las=2)
axisPhylo(side=1)
dev.off()


#par(mfrow=c(2,1))
#par(mar=c(4,4,4,4))
#plotBS(midpoint(fitGTR$tree), bs.GTR, p = 90, type="p")
#title("a)")
#cnet <- consensusNet(bs.GTR, p=0.2)
#plot(cnet, "2D", show.edge.label=TRUE)
#title("b)")

# comparing a phylogeny without evol model,
# with one produced the diff transitions
comparePhylo(small_NJ, fitGTRopt.RodMar$tree, 
             plot = TRUE, force.rooted = TRUE)

## depth of correction
dist_mat <- as.matrix(dist.p(gen_small, 
                             cost = "polymorphism", 
                             ignore.indels = TRUE),diag=T)

## mantel test testing p-distance and the distance generated by the GTR model 

mantel(dist_mat,
       as.dist(cophenetic(fitGTRopt.RodMar$tree),diag=T))

## checking correspondence between distances
x <- as.vector(as.dist(cophenetic(small_NJ)))
y <- as.vector(as.dist(cophenetic(fitGTRopt.RodMar$tree)))
plot(x, y, xlab="Pairwise distances on NJ tree", 
     ylab="Pairwise distances on the optimeed tree (GTR model)", 
     main="Is NJ appropriate?", 
     pch=20,  cex=1,col="gray")
abline(lm(y~x), col="red",lwd=2)
cor(y,x)^2



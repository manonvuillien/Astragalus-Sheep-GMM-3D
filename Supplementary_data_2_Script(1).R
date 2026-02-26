################################################################################
############# SCRIPT Exploring biological and ecological components ############
############# of sheep astragalus size and shape variation using 3D ############
############# geometric morphometrics: towards a bioarchaeological proxy #######
########### VUILLIEN et al., Journal of Archaeological Method and Theory #######
################################################################################

################################################################################
################ IMPORTATION DES LIBRAIRIES NECESSAIRES AU CODE ################
################################################################################

library(devtools)
library(stats)
library(tidyr)
library(geomorph)
library(Morpho)
library(rgl)

################################################################################
##################### FONCTIONS UTILES POUR L'AUTOMATISATION ###################
################################################################################

# function resampling a single curve
cursub.interpo <- function(cur,req){
  mat = as.matrix(dist(cur))
  DPO = NULL
  for(j in 1:nrow(cur)-1){
    a = mat[j+1,j]
    DPO = c(DPO,a)}
  
  DCO = NULL
  for(j in 1:length(DPO)){DCO[j]=sum(DPO[1:j])}
  DN = sum(DPO) / (req+1)
  DCN = DN*(1:(req))
  proxima = matrix(nrow=req,ncol=2)
  
  for(k in 1:length(DCN)){
    first = which.min(abs(DCO- DCN[k]))
    proxima[k,1] = first
    second = which.min(abs(DCO[-first]- DCN[k]))
    ifelse(first==second,second<-(second+1),second<-second)
    proxima[k,2] = second}
  proxima = proxima + 1 
  proxima2 = matrix(nrow=req,ncol=2)
  
  for(i in 1:req){
    if(proxima[i,1] > proxima[i,2]){
      proxima2[i,1]=proxima[i,2];proxima2[i,2]=proxima[i,1]}
    else if(proxima[i,1] < proxima[i,2]){
      proxima2[i,1]=proxima[i,1];proxima2[i,2]=proxima[i,2]}}
  
  VEC = matrix(nrow = req, ncol = 3)
  for(l in 1:req){VEC[l,]=as.matrix(cur[proxima2[l,2],]-cur[proxima2[l,1],])}
  COMP = NULL
  for(n in 1:req){COMP[n]=DCN[n]-DCO[proxima2[n,1]-1]}
  NOR = NULL
  for(m in 1:req){NOR[m]=COMP[m]/(DPO[proxima2[m,1]])}
  
  VECF = VEC * NOR
  PTS = cur[proxima2[,1],] + VECF
  PTS = rbind(cur[1,], PTS, cur[dim(cur)[1],])
  PTS
}

subsampl.inter <- function(matlm,curlist,required,fix){
  if(is.list(curlist)==F)
    print("curlist must be a list giving the curves(rowindex)")
  else 
    if(is.vector(required)==F)
      print("required must be a vector giving the number of points required per curve")
  else 
    if(length(curlist)!=length(required))
      print("curlist and required must be of same length")
  else 
    output=matlm[fix,]
  for (i in 1:length(curlist)){
    cur = matlm[curlist[[i]],]
    req = required[i]
    out = cursub.interpo(cur,req)
    rownames(out) = paste("curve", i, "-", (1:dim(out)[1]-1), sep="")
    output = rbind(output, out[2:(dim(out)[1]-1),])
  }
  output
}

################################################################################
######################### IMPORTATION DU JEU DE DONNEES ########################
################################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

talus = read.csv("Template_Evosheep_subset.csv", header=TRUE, sep=";", dec=".")
talus_data = talus[,2:dim(talus)[2]]
talus_names = talus$Specimen.mm.

################################################################################
################################ VARIABLES FIXES ###############################
################################################################################

nbr_talus = dim(talus)[1] - 1
n = nrow(talus_data) # Nombre d'individus
tem = t(talus_data)
dimx = 3 # nombre de colonnes
npts = dim(talus_data)[2]/dimx # nombre de points par colonne

################################################################################
################################## APPLICATION #################################
################################################################################

# Calcul
tem = array(tem, c(dimx, npts, n))
tem = aperm(tem, c(2, 1, 3))

curlist = list(c(12:84), c(85:121))
required = c(length(curlist[[1]]),length(curlist[[2]])) 
fix = (1:11)
req = 198
mesh.atlas = read.ply("LBOAWK001_talus.ply", ShowSpecimen = FALSE, addNormals = TRUE)

list_curves = c()
list_atlas = c()
list_mesh = c()
data = c()
list_names_mesh = talus_names

# Rééchantillonnage de toutes les courbes
list_curves = c()
for(i in 1:5){list_curves[[i]] = subsampl.inter(tem[,,i], curlist, required, fix)}

# Création de l'individu de référence
patch = tem[,,1][122:198,]
indref = rbind(list_curves[[1]], patch)
rownames(indref) = c(1:198)
  
# Création de l'atlas de référence
atlas = createAtlas(mesh.atlas, landmarks = indref[c(1:11),], patch = indref[-c(1:11),], patchCurves = list(c(12:84),c(85:121)))
plotAtlas(atlas)
mesh2ply(mesh.atlas, filename=list_names_mesh[[1]])

data[[1]] = indref
##keep.fix = c(1:11)

# Création du mesh de référence
list_mesh[[1]] = mesh.atlas

for(j in 1:nbr_talus+1){
  print(paste0("Mesh n° ",j))
  
  patch = tem[,,j][122:198,]
  ind_cible = rbind(list_curves[[j]], patch)
  rownames(ind_cible) = c(1:198)
  
  mesh_cible = tps3d(mesh.atlas,indref,ind_cible,threads=1)
  list_mesh[[j]] = mesh_cible
  #plot3d(mesh_cible,type="wire", size=1, add=T, col="#99CCFF")
    
  # Enregistrement des meshs
  mesh2ply(mesh_cible, filename = list_names_mesh[[j]])
  data[[j]] = ind_cible
  }

surp = c(1:nrow(indref))[-fix]
data = bindArr(data, along = 3)
dimnames(data)[[3]] = list_names_mesh

slide = slider3d(data, SMvector = fix, deselect = TRUE, surp = surp, iterations = 1, 
                 meshlist = list_mesh, mc.cores = 1, fixRepro = TRUE, sur.path = ".")

dataslide = c()
for(i in 0:nbr_talus+1){
  dt = data.frame(data.frame(slide$dataslide[,,i]) %>% pivot_longer(cols=c(1,2,3),values_to='points'))
  dt = t(dt$points)
  dataslide = data.frame(rbind(dataslide, c(list_names_mesh[i] , dt)))}
colnames(dataslide) = c("ID",paste0(rep(c("X", "Y", "Z"), length.out = (ncol(dataslide)-1)), rep(1:((ncol(dataslide)-1)/3), each=3)))
write.csv(dataslide, "data_slide_sheeps_talus.csv")

################################################################################
############# END - SCRIPT Exploring biological and ecological components ############
############# of sheep astragalus size and shape variation using 3D ############
############# geometric morphometrics: towards a bioarchaeological proxy #######
########### VUILLIEN et al., Journal of Archaeological Method and Theory #######
################################################################################
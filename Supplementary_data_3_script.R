################################################################################
############# SCRIPT Exploring biological and ecological components ############
############# of sheep astragalus size and shape variation using 3D ############
############# geometric morphometrics: towards a bioarchaeological proxy #######
########### VUILLIEN et al., Journal of Archaeological Method and Theory #######
################################################################################

################################################################################
################################ LIBRARIES #####################################
################################################################################

library(geomorph) 
library(Morpho) 
library(corrplot) 
library(ggplot2) 
library(rstatix) # to open the dataset
library(rgl) # window for 3D visualisation

################################################################################
################################# DATA #########################################
################################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

evo_talus = data.frame(read_delim("Supplementary_data_3_data_sheep.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE), check.names = FALSE)
evo_talus[, 16:610] = lapply(evo_talus[, 16:610], as.numeric)
evo_talus$Mobility = as.character(evo_talus$Mobility)

################################################################################
############################# GROUPS  ##########################################
################################################################################

gp1 = evo_talus$Breed
gp2 = evo_talus$Species
gp4 = evo_talus$Sexe[-c(74:75, 88:96)]
gp5 = evo_talus$Weight[c(1:18, 42:52)]
gp6 = evo_talus$Length[c(1:18, 22:25, 34:37)]
gp7 = evo_talus$Horn[-c(61, 74:75, 56, 88:96)]
gp8 = evo_talus$Coat[-c(74:96)]
gp9 = evo_talus$Tail[-c(74:96)]
gp10 = evo_talus$Country[-c(74:75, 88:96)]
gp11 = evo_talus$`Classification KG`[-c(74:96)]
gp12 = evo_talus$`Climatic zone`[-c(74:75, 88:96)]
gp13 = evo_talus$Elevation[-c(74:96)]
gp14 = evo_talus$Mobility[-c(74:96)]

################################################################################
####################################### PCA ####################################
################################################################################

pcoord = arrayspecs(evo_talus[,17:610], p = 198, k = 3) 
sup = gpagen(pcoord, PrinAxes = TRUE) # Procruste superimposition
plotOutliers(sup$coords, inspect.outliers = TRUE)

PCA = gm.prcomp(sup$coords)
dt = summary(PCA)[[1]]["Cumulative Proportion",]
plot(1:95,round(dt,2)) # the number of components depending on the dataset 
abline(v = which(round(dt,2) >= 0.8)[1], h = 0.8)

newdata = data.frame(cbind(evo_talus[,-c(17:610)], PCA$x[,1:14]), check.names = FALSE)

##Plot PCA
couleurs = c("darksalmon", "darkorchid", "darkred", "aquamarine3","blue","cornflowerblue", "grey", "forestgreen","green", "chartreuse", "deepskyblue", "black","brown", "darkturquoise", "darkorange")
base = ggplot(newdata, aes(Comp1, Comp2, colour = `Breed`)) +
  geom_point(size = 3.5, aes(shape = `Breed`)) +
  xlab(paste0("PC1 - ", round(summary(PCA)[[1]]["Proportion of Variance", 1] * 100, 2), "%")) +
  ylab(paste0("PC2 - ", round(summary(PCA)[[1]]["Proportion of Variance", 2] * 100, 2), "%")) +
  scale_color_manual(values = couleurs) +
  scale_shape_manual(values = c(16, 16, 15, 17, 17, 17, 18, 18, 18, 18, 17, 3, 4, 17, 15)) +
  theme_bw() + 
  theme(legend.background = element_rect(fill = "white", linewidth = 0.5, colour = "grey70", linetype = "solid"), 
        axis.text = element_text(size = 10, face = "bold"), 
        aspect.ratio = 1, 
        panel.grid = element_blank())

hull = data.frame(newdata %>% group_by(Breed) %>% slice(chull(Comp1, Comp2)), check.names = FALSE)
base + geom_polygon(data = hull, aes(x = Comp1, y = Comp2, fill = Breed, color = Breed), alpha = 0)

##Plot PCA with centroïds - Mandak and Kuyuruksuz were excluded for the visualisation (Figure 5 in the present article)
couleurs = c("darksalmon", "darkorchid", "darkred", "aquamarine3","blue","cornflowerblue", "grey", "forestgreen", "deepskyblue", "black","brown", "darkturquoise", "darkorange")
newdata_filtered <- newdata %>%
  filter(!Breed %in% c("Mandak", "Kuyuruksuz"))

centroids <- newdata_filtered %>%
  group_by(Breed) %>%
  summarise(
    c1 = mean(Comp1),
    c2 = mean(Comp2),
    .groups = "drop"
  )
hull <- newdata_filtered %>%
  group_by(Breed) %>%
  slice(chull(Comp1, Comp2))

ggplot(newdata_filtered, aes(Comp1, Comp2, colour = Breed)) +
  
  geom_polygon(data = hull, aes(fill = Breed), alpha = 0.08, colour = NA) +
  geom_point(size = 3, alpha = 0.6) +
  
  geom_point(data = centroids,
             aes(x = c1, y = c2, fill = Breed),
             shape = 21, size = 7, colour = "black", stroke = 1.5, inherit.aes = FALSE) +
  
  scale_color_manual(values = couleurs) +
  scale_fill_manual(values = couleurs) +
  
  xlab(paste0("PC1 - ", round(summary(PCA)[[1]]["Proportion of Variance", 1] * 100, 2), "%")) +
  ylab(paste0("PC2 - ", round(summary(PCA)[[1]]["Proportion of Variance", 2] * 100, 2), "%")) +
  
  theme_bw() +
  theme(
    legend.background = element_rect(fill = "white", linewidth = 0.5, colour = "grey70"),
    axis.text = element_text(size = 10, face = "bold"),
    aspect.ratio = 1,
    panel.grid = element_blank()
  )


################################################################################
###################################### SIZE ####################################
################################################################################

size = sup$Csize

newdataCS = cbind(evo_talus[,-c(17:610)], log(size)) #dataset 
shapiro.test(newdataCS$`log(size)`) # normality test
ggqqplot(log(size)) # normality visualisation
bartlett.test(`log(size)`~Species, data=newdataCS) # homogeneity of variances

# Change groups for testing other qualitative variables
group = gp1[-c(74:87)] #remove Mandak, Kuyuruksuz, Turkish breed unidentified specimens
kruskal.test(log(size[-c(74:87)]) ~ group) # group 
CS = pairwise.wilcox.test(log(size), group, p.ajust = "bonferroni", na.rm = TRUE)
write.csv(CS$p.value, file = "CS_pairwise_Wilcoxon_Breed.csv")

#Group changes: 
kruskal.test(log(size[-c(61, 74:75, 56, 88:96)]) ~ gp7) # test anova Horn

# Testing relation between Size, Sexe and Presence or absence of horn
anova_test(`log(size)` ~Sexe+Horn, data = na.omit(newdataCS[,c("log(size)", "Sexe", "Horn")]))

# Plot CS Size ~ Horn VS Sexe
dt = na.omit(newdataCS[,c("Horn","log(size)","Sexe")])
dt = dt[dt$Horn %in% c("No","Yes"), ]
ggplot(dt, aes(Horn, `log(size)`, fill = Sexe)) +
  geom_boxplot(aes(group = interaction(Horn, Sexe)), width = 0.5, alpha = 1, position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(fill = Sexe), color = "black", shape = 21, position = position_dodge(width = 0.75), size = 2.5, alpha = 0.8, stroke = 0.6) + 
  scale_color_manual(values = c("#A6CEE3", "#1F78B4")) + 
  scale_fill_brewer(palette = "Paired") + theme_bw() + 
  theme(legend.background = element_rect(fill = "white", linewidth = 1, colour = "grey70", linetype = "solid"),
        axis.text = element_text(size = 12), panel.grid = element_blank())

# Plot CS Size ~ Sheep Breeds/landraces VS Country
newdataCS$Breed = factor(newdataCS$Breed, levels = c("Awassi_AB","Awassi_WK","Bakhtiari_LL","Ziaran","Inconnue","Karaksz",
                                                     "Kuyuruksuz","Mandak","Gumz","Bonga","Washera","Menz","Black Headed Somali",
                                                     "Ovis_gmelini_musimon","Ovis_gmelini"))
newdataCS_Country = newdataCS[,c("Breed","log(size)","Country")]

couleurs = c("darksalmon","darkorchid","darkred","darkorange", "grey", "forestgreen","green", "chartreuse","cornflowerblue","blue","darkturquoise", "deepskyblue","aquamarine3", "black","brown")
ggplot(na.omit(newdataCS_Country[,c("Breed","log(size)","Country")]), 
       aes(x = Breed, y = `log(size)`, fill = factor(Breed))) +
  geom_boxplot(aes(group = interaction(Breed, Country)), width = 0.5, alpha = 1, position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(color = Country), position = position_dodge(width = 0.75), size =2.5, alpha = 0.8) +
  scale_fill_manual(values = couleurs) + scale_color_manual(values = rep("black", length(unique(newdataCS_Country$Country)))) + 
  guides(fill = "none") + theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 1, colour = "black", linetype = "solid"),
        axis.text = element_text(size = 12, angle = 40), panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 40, vjust = 1, hjust = 1))+
  labs(title = "Size by Breed and Country", x = "Breed", y = "log(size)")

# Plot CS Size ~ Country # Supplementary data 4
couleurs = c("darksalmon", "darkorchid", "darkred", "aquamarine3","blue","cornflowerblue", "grey", "forestgreen","green", "chartreuse", "deepskyblue", "black","brown", "darkturquoise", "darkorange")
ggplot(na.omit(newdataCS[,c("Country","log(size)")]), aes(Country, `log(size)`, fill=Country)) + geom_boxplot(position = position_dodge(1), aes(shape=Country)) +
  geom_dotplot(binaxis ='y', stackdir = 'center', binwidth=0.01,  dotsize = 0.75, binpositions ="all", color="black") +
  scale_fill_manual(values = couleurs) + theme_bw() + 
  theme(legend.background = element_rect(fill = "white", linewidth = 1, colour = "grey70", linetype = "solid"),
        axis.text = element_text(size=12), panel.grid = element_blank())

################################################################################
################################### CORRELATIONS ###############################
################################################################################

mycor <- function(x,cor_method){
  r = apply(x, 2, function(j){apply(x, 2, function(i){as.numeric(cor.test(i,j,method = cor_method)$estimate)})})
  P = apply(x, 2, function(j){apply(x, 2, function(i){as.numeric(cor.test(i,j,method = cor_method)$p.value)})})
  out = c()
  out$P = P
  out$r = r
  return(out)}

quant = sapply(newdataCS, function(x) is.numeric(x) | is.integer(x))
data_corr = newdataCS[ ,quant]
corr = mycor(data_corr, "spearman")
corr

corrplot(corr = corr$r, mar = c(0,0,1,0), method = "shade", type = "lower", sig.level = 0.01, addCoef.col = "black", 
         insig = "pch", pch.cex = 0, pch.col = "black", tl.col = "black",col = COL2("RdBu",200))


################################################################################
##################################### MANOVA ###################################
################################################################################

# Application example of 3.3 "Effect of breeding selection on sheep astragalus size and shape variation"

# For Breed/Species ANOVA
group = gp1

# The categories removed from the following analyses are: "Mandak", "Kuyuruksuz", "Inconnue", "Ovis_gmelini"
# They correspond to the following indices: -c(74:87, 95:96)

group = group[-c(74:87, 95:96)]
sup2 = sup$coords[,,-c(74:87, 95:96)]
supsize = sup$Csize[-c(74:87, 95:96)]

gdf = geomorph.data.frame(shape = sup2, size = supsize, taxa = group)
fit = procD.lm(shape ~taxa, data = gdf, iter = 999, RRPP = TRUE, na.omit = T) # RRPP = TRUE randomize residuals
fit$aov.table
summary(fit$residuals)

# PC plot rotated to major axis of fitted values
plot(fit, type = "diagnostics", pch = 19, col = "blue", reg.type = c("RegScore"), outliers = FALSE) 
# Extracting objects and plotting options
# diagnostic plots, including plotOutliers
plot(fit, type = "diagnostics", outliers = TRUE) 

PW = pairwise(fit, groups = group, covariate = NULL)
summary(PW)
summary(fit$residuals)

# For all other variables, change here according to their information group
group2 = evo_talus$`Classification KG` #example here for KG classification group
group2 = group2[-c(74:96)]
sup2 = sup$coords[,,-c(74:96)]
supsize = sup$Csize[-c(74:96)]

gdf = geomorph.data.frame(shape = sup2, size = supsize, taxa = group2)
fit = procD.lm(shape ~taxa, data = gdf, iter = 999, RRPP = TRUE, na.omit = T) # RRPP = TRUE randomize residuals
fit$aov.table
summary(fit$residuals)

# PC plot rotated to major axis of fitted values
plot(fit, type = "diagnostics", pch = 19, col = "blue", reg.type = c("RegScore"), outliers = FALSE) 
# Extracting objects and plotting options
# diagnostic plots, including plotOutliers
plot(fit, type = "diagnostics", outliers = TRUE) 

PW = pairwise(fit, groups = group2, covariate = NULL)
summary(PW)
summary(fit$residuals)

################################################################################
#################### 3D PROJECTION only for sheep dataset ######################
################################################################################

mesh = read.ply("Supplementary_data_3_Ziaran_ovis_MMM_test.ply", ShowSpecimen = FALSE, addNormals = TRUE)
meshcoord = sup$coords[,,findMeanSpec(sup$coords)]

PCA_1_mesh = mshape(PCA$shapes$shapes.comp1$min)
PCA_2_mesh = mshape(PCA$shapes$shapes.comp1$max)
msh1 = tps3d(mesh, meshcoord, PCA_1_mesh)
msh2 = tps3d(mesh, meshcoord, PCA_2_mesh)

require(rgl)
mfrow3d(1, 2, sharedMouse = T)
shade3d(x = msh1, col = "grey70", lit = T)
next3d()
shade3d(x = msh2, col = "#99CCFF", lit = T)
rglwidget()

################################################################################
################################################################################
############# SCRIPT Exploring biological and ecological components ############
############# of sheep astragalus size and shape variation using 3D ############
############# geometric morphometrics: towards a bioarchaeological proxy #######
########### VUILLIEN et al., Journal of Archaeological Method and Theory #######
################################################################################
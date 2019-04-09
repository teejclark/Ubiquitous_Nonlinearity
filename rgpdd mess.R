#### Animal Population Dynamics Code
#### T.J. Clark - 5/20/18

rm(list=ls(all=T))

library("rgpdd")
library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")

###################################################################

#### DATA PREPARATION

setwd("~/PhD/R Code/Animal Population Dynamics")
extra <- read.csv("Extra Data.csv")
extra <- extra[,1:20]
extra <- extra[1:4769,] # remember to revise later

#### Get GPDD Data; Remove extraneous columns
gpdd_join <- gpdd_main %>%
  inner_join(gpdd_taxon,by="TaxonID") %>%
  inner_join(gpdd_data,by="MainID") %>%
  select(-DataSourceID, -TaxonID,-BiotopeID, -LocationID,
         -SourceDimension, -SpatialDensity, #-SourceTransform,
         -SourceTransformReference,-AssociatedDataSets,
         -SamplingFrequency, -SiblyFittedTheta,
         -SiblyThetaCILower, -SiblyThetaCIUpper,
         -SiblyExtremeNEffect, -SiblyReturnRate,
         -SiblyCarryingCapacity, -WoldaCode,
         -Authority, - Notes.y, - DataID, -PopulationUntransformed,
         -TimePeriodID, - Generation, - SeriesStep,
         -DecimalYearBegin, - DecimalYearEnd)

# Join with Gray Data
gpdd_join <- rbind(gpdd_join,extra)

# Join to the Life History Dataset
lifehistory <- read.csv("Life History Data.csv")

lifehistory <- lifehistory %>%
  unite("TaxonName", c("Genus","Species"),
        sep=" ",remove=T)

# join together stuff with life history data, leave out missing
gpdd_join <- gpdd_join %>%
  inner_join(lifehistory, by="TaxonName")


#### FILTERING DATA

# STEP 1: Remove t<30
# STEP 2: Remove "unreliable databases" <3
# STEP 3: Include only bonyfish,insects,mammals,birds...
# STEP 4: Filter out negative values and all 0s
# STEP 5: Remove those with less than 5 unique values
# STEP 6: (Remove those with more than XX repeating zeros)


# function to remove 0s
zero_fxn <- function(x){
  k<-rle(x$Population)
  m<-k$lengths[k$values==0]
  ifelse(length(max(m))>0,
         max(m),0)}

# Filter
gpdd_filter <- gpdd_join %>%
  filter(Reliability>2)%>% # remove unreliable datasets
  filter(TaxonomicClass!="Bacillariophyceae") %>% #remove other animals
  filter(TaxonomicClass!="Dinophyceae")%>%
  filter(TaxonomicClass!="Bivalvia") %>%
  filter(TaxonomicClass!="Chondrichthyes") %>%
  filter(TaxonomicClass!="Crustacea") %>%
  filter(TaxonomicClass!="Echinoidea") %>%
  filter(TaxonomicClass!="Polychaeta") %>%
  filter(TaxonomicClass!="Scyphozoa") %>%
  filter(TaxonomicClass!="Gastropoda") %>%
  filter(TaxonomicClass!="Unknown") %>%
  filter(TaxonomicClass!="Hyperoartia") %>%
  group_by(MainID) %>%
  mutate(Population=replace(Population,SourceTransform=="Log",
                            10^Population)) %>% #back-transform log
  mutate(a = sum(Population), # can't be all zeros
         b=min(Population), # no negative abundances
         c=n(), # datset length => 30
         d=n_distinct(Population), # unique values
         e=zero_fxn(.data)) %>% # remove repeating zeros
  filter(a>0) %>% filter(b>=0) %>% filter(c>=30) %>%
  filter(d>=5) %>%
  #filter(e<=15|is.na(e)) %>%
  ungroup()

# gpdd_filter %>% group_by(MainID) %>%
#   filter(TaxonomicClass=="Insecta") %>%
#   ggplot(aes(x=e))+geom_bar()
#
# gpdd_filter %>% group_by(TaxonomicClass) %>%
#   summarize(Length=length(unique(MainID)),
#             Species=length(unique(TaxonName)))


# Attach Corrected Data - IMPORTANT
correct <- read.csv("Missing Data_big.csv")
# stick corrected life history data onto dataset
gpdd_filter <- rbind(gpdd_filter,correct)

#### FIRST-DIFFERENCE//SCALE
gpdd_st <- gpdd_filter %>%
  group_by(MainID) %>%
  mutate(first_diff=Population-lag(Population)) %>% #1st diff
  slice(-1) %>%
  mutate(sc_Population=scale(first_diff)) %>% ungroup() #scale


###############################################################
#### BOX 1

## Extra: anti-join to find missing life history data
## attach to good dataset with inner join
anti_filter <- gpdd_filter %>%
  anti_join(lifehistory,by="TaxonName")

write.csv(anti_filter,"Missing Data_big.csv")
###############################################################

#### FUNCTIONS

# FUNCTION 1 = CALCULATE E - DIMENSIONALITY
simplex_fun<-function(x){
  which.max(unlist(
    simplex(time_series = x$sc_Population,
                  lib= c(1,length(x$sc_Population)),
                  pred= c(1,length(x$sc_Population)),
                  tau=1)["rho"]))
}

# EXTRA FUNCTION 1.5 - Calculate E and Tau
simplex_extra_fun<-function(x){
  output <-   simplex(time_series = x$sc_Population,
            lib= c(1,length(x$sc_Population)),
            pred= c(1,length(x$sc_Population)),
            E=1:10,tau=1:10)
  output$E[which.max(output$rho)]
}


# FUNCTION 2 = CALCULATE THETA - NONLINEARITY

theta_fun <- function(x){
 output <- s_map(time_series = x$sc_Population,
        norm_type = "L2 norm",
        lib =  c(1,length(x$sc_Population)),
        pred = c(1,length(x$sc_Population)),
        E = min(x$E))

output$theta[which.max(output$rho)]
}

# FUNCTION 2.5 = CALCULATE NONLINEARITY FROM THETA's MAE
# randomization procedure to categorize nonlinearity

nonlin_fun <- function(x){
  output <- s_map(time_series = x$sc_Population,
                  norm_type = "L2 norm",
                  lib =  c(1,length(x$sc_Population)),
                  pred = c(1,length(x$sc_Population)),
                  E = min(x$E))

  output$mae[which(output$theta==0)]-min(output$mae)
}


# FUNCTION 3 = CALCULATE Tp - FORECAST SKILL

forecast_fun <- function(x) {
  output <- s_map(time_series = x$sc_Population,
        norm_type = "L2 norm",
        lib = c(1,length(x$sc_Population)),
        pred = c(1,length(x$sc_Population)),
        E = min(x$E),
        theta = min(x$theta),
        tp = c(seq(1,5,1)))

  output$rho[which(output$tp==1)]
}

#FUNCTION 4 = Calculate p values for forecast skill
pred_fun <- function(x) {
  output <- s_map(time_series = x$sc_Population,
                  norm_type = "L2 norm",
                  lib = c(1,length(x$sc_Population)),
                  pred = c(1,length(x$sc_Population)),
                  E = min(x$E),
                  theta = min(x$theta),
                  tp = c(seq(1,5,1)))

  output$p_val[which(output$tp==1)]
}

###################################################################

#### MODEL RUNS

# CALCULATE E and Theta; also CV, and Time-mismatch
# to calculate CVs to not have -INF; change 0 to 0.001
gpdd_1 <- gpdd_st %>%
  group_by(MainID) %>%
  mutate(first_diff_zeroes=replace(first_diff,first_diff==0,0.001)) %>%
  mutate(CV=sd(first_diff_zeroes/mean(first_diff_zeroes))) %>%
  mutate(E=simplex_extra_fun(.data)) %>%
  mutate(theta=theta_fun(.data)) %>%
  mutate(time_mismatch=c/MinAge)

#### RANDOMIZATION TEST - non-linear for-loop

# original function - to calculate original mae's
gpdd_2 <- gpdd_1 %>%
  group_by(MainID) %>%
  summarize(mae=nonlin_fun(.data))

#####################################################################
#### Box 2
# randomization test - calculate other mae's
reps <- 1000
mat <- matrix(0,reps,length(unique(gpdd_1$MainID)),
              dimnames=list(seq(1,reps,1),unique(gpdd_1$MainID)))

for (i in 1:reps){

  x <-  gpdd_1 %>%
    group_by(MainID) %>%
  .[sample(1:nrow(.)),] %>%
   summarize(mae=nonlin_fun(.data))

  mat[i,] <- x$mae

}
#####################################################################

#load randomization dataset
random <- read.csv("random.csv")

# get the 95% quantiles of all randomization simulations
conf <- data.frame(lapply(random,quantile,probs=0.95,na.rm=T))
c <- data.frame(t(conf))

# summarize original MAEs
gpdd_2.5 <- gpdd_2 %>% group_by(MainID) %>%
  summarize(mae=min(mae))
# combine dataframe - original and randomized MAEs
blah <- data.frame(gpdd_2.5,c)

# if statement to label mae's
# if original mae > randomization = nonlinear
gpdd_2.5_fin <- blah %>%
  mutate(nonlin=ifelse(
    mae>X95.,
    "nonlinear",
    "linear"))


# Calculate predictability and p-value
gpdd_3 <- gpdd_1 %>%
  group_by(MainID) %>%
  mutate(rho=forecast_fun(.data)) %>%
  mutate(p=pred_fun(.data))

# combine all final stuff
# also change sampling unit names
gpdd_fin <- left_join(gpdd_3,gpdd_2.5_fin) %>%
  mutate(TaxonomicClass=replace(TaxonomicClass,
                                TaxonomicClass=="Actinopterygii",
                                "Osteichthyes")) %>%
  droplevels()

# a little editing for sampling stuff
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Adults"] <- "Breeding individuals"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Breeding pairs"] <- "Breeding individuals"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Eggs"] <- "Young"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Females"] <- "Breeding individuals"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Host plant damaged"] <- "Individuals killed"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Individuals caught"] <- "Trapped Individuals"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Individuals exported"] <- "Individuals killed"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Kg (females)"] <- "Biomass"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Male territories or breeding pairs"] <- "Breeding individuals"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Males"] <- "Breeding individuals"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Nests"] <- "Breeding individuals"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Pupae"] <- "Young"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Pairs"] <- "Breeding individuals"
gpdd_fin$SamplingUnits[gpdd_fin$SamplingUnits=="Breeding Pairs"] <- "Breeding individuals"


#####################################################################
#### BOX 3
# test function
test_fun<-function(x){

  ifelse(length(which.max(unlist((simplex(ts <- x$Population,
                                          lib<- c(1,length(x$Population)),
                                          pred <- c(1,length(x$Population)),
                                          tau=1)["rho"]))))==0,
         stop(print(x$MainID)),
         which.max(unlist((simplex(ts <- x$Population,
                                   lib<- c(1,length(x$Population)),
                                   pred <- c(1,length(x$Population)),
                                   tau=1)["rho"]))))
}
#####################################################################

#####################################################################
#### SUMMARY STATS

# TOTAL
# 705 datasets; 632 predictable // 201 species; 159 predictable
# BREAKDOWN BY CLASS - using t>=30
# Birds: 70 datasets; 38 predictable // 53 species; 32 predictable
# Insects: 523 datasets; 504 predictable // 46 species; 42 predictable
# Mammals: 35 datasets; 24 predictable// 25 species; 19 predictable
# Bony Fish: 77 datasets; 66 predictable// 77 species; 66 predictable

gpdd_fin %>% group_by(TaxonomicClass) %>%
  summarize(Length=length(unique(MainID)),
            Species=length(unique(TaxonName)))


# PCA STUFF - remember - use only for analysis -
# because I'm dropping NA datsets
# for visualization
library(factoextra)
library(corrplot)

gpdd_fin <- gpdd_fin_yaya

# remove nas
gpdd_fin_blah <- gpdd_fin[!is.na(gpdd_fin$Len),]
gpdd_fin_blah <- gpdd_fin_blah[!is.na(gpdd_fin_blah$MinAge),]
gpdd_fin_blah <- gpdd_fin_blah[!is.na(gpdd_fin_blah$Lifesp),]
gpdd_fin_blah <- gpdd_fin_blah[!is.na(gpdd_fin_blah$Fert),]

# repeat for p
gpdd_p_blah <- gpdd_fin %>%
  filter(p<= 0.05)
gpdd_p_blah <- gpdd_p_blah[!is.na(gpdd_p_blah$Len),]
gpdd_p_blah <- gpdd_p_blah[!is.na(gpdd_p_blah$MinAge),]
gpdd_p_blah <- gpdd_p_blah[!is.na(gpdd_p_blah$Lifesp),]
gpdd_p_blah <- gpdd_p_blah[!is.na(gpdd_p_blah$Fert),]


#####################################################################
#### Here's another FREAKING Box
#### Phylogenetic Correction before PCA - Revell (2009)
library(ape)
library(phytools)

a_old <- read.tree("newick.nwk") # old tree
a <- read.tree("animal_tree.nwk") # better tree

# load datset and shape up
anim_data <- read.csv("anim_data_fin.csv")
anim_data <- anim_data[-152,]
anim_data <- anim_data[-152,]
anim_data <- anim_data[-217,]
anim_data <- anim_data[-126,]
anim_data <- anim_data[-127,]

anim_data <- data.frame(anim_data,row.names=anim_data$CorrectName)

# snip tree
a <- drop.tip(a,c("Callionymus_lyra","Liparis_liparis",
                  "Liza_ramado","Mullus_surmuletus"))

# rename
gpdd_fin_core <- anim_data

# create PCA matrix
X <- as.matrix(data.frame(log(gpdd_fin_core$Len),
                          log(gpdd_fin_core$MinAge),
               log(gpdd_fin_core$Lifesp),log(gpdd_fin_core$Fert)))
rownames(X) <- rownames(gpdd_fin_core)

# function to find phylogenetic eigenvalues
# gonna run manually to troubleshoot
pca_results <- phyl.pca(a,X)
biplot(pca_results)

gpdd_fin_core <- data.frame(cbind(gpdd_fin_core,pca_results$S))

#####################################################################

# store data and run PCA
res.pca <- prcomp(cbind(log(gpdd_fin_blah$Len),
                                log(gpdd_fin_blah$MinAge),
                                log(gpdd_fin_blah$Lifesp),
                                log(gpdd_fin_blah$Fert)),scale=T)
res.pca.p <- prcomp(cbind(log(gpdd_p_blah$Len),
                        log(gpdd_p_blah$MinAge),
                        log(gpdd_p_blah$Lifesp),
                        log(gpdd_p_blah$Fert)),scale=T)

# visualize
windows()
fviz_pca(res.pca.p,
         geom.ind="point",
         col.ind=gpdd_p_blah$TaxonomicClass,
         palette=c("#F0E442","#0072B2", "#D55E00", "#CC79A7"),
         addEllipses=T,
         legend.title="Groups",xlab="PC1",ylab="PC2")+
  theme(text=element_text(size=20))

# cbind to original data
res.pca.blah <- data.frame(res.pca$x)
gpdd_fin <- bind_cols(gpdd_fin_blah,res.pca.blah)

res.pca.p.blah <- data.frame(res.pca.p$x)
gpdd_p <- bind_cols(gpdd_p_blah,res.pca.p.blah)

# Remove redundant time-series
gpdd_fin_core <- gpdd_fin %>%
  group_by(TaxonName) %>%
  sample_n(1)

write.csv(gpdd_fin_core,"anim_data_fin.csv")


#### PHYLOGENETIC STUFF
library(ape)
library(adephylo)
library(nlme)
library(geiger)

#####################################################################
#### BOX 77 - Creating Phylogenetic Tree - ONE VARIABLE
#### this is garbage remove eventually

a <- read.tree("newick.nwk")
b <- as.matrix(distTips(a,tips="all",method="Abouheif"))
c <- data.frame(b[6,])
d <- rownames_to_column(c)
write.csv(d,"phylo.data.csv")

# okay now add in values
phylo_data <- read.csv("phylo.data.csv")
# combine
gpdd_fin_core <- bind_cols(gpdd_fin_core,phylo_data) %>%
  select(-X, -rowname, -b.6..., -X.1,-old.name)

# get p values
gpdd_p_core <- gpdd_fin_core %>%
  filter(p<=0.05)

# also calculate p
# gpdd_p_core <- gpdd_p %>%
#   group_by(TaxonName) %>%
#   sample_n(1)

#####################################################################

#### Create Phylogenetic Variance-Covariance matrix
a <- read.tree("newick.nwk")# old tree
a <- read.tree("animal_tree.nwk")

# load datset and shape up
anim_data <- read.csv("anim_data_fin.csv")
anim_data <- anim_data[-152,]
anim_data <- anim_data[-152,]
anim_data <- anim_data[-217,]
anim_data <- anim_data[-126,]
anim_data <- anim_data[-127,]

anim_data <- data.frame(anim_data,row.names=anim_data$CorrectName)

# snip tree
a <- drop.tip(a,c("Callionymus_lyra","Liparis_liparis",
                  "Liza_ramado","Mullus_surmuletus"))

# check that data matches
obj <- name.check(a,anim_data)

# rename
gpdd_fin_core <- anim_data

# p
gpdd_p_core <- gpdd_fin_core %>%
  filter(p<=0.05)
gpdd_p_core <- data.frame(gpdd_p_core,row.names=gpdd_p_core$CorrectName)

# remove tips that don't match with data
obj_2 <- name.check(a,gpdd_p_core)
b <- drop.tip(a,c("Acinonyx_jubatus","Alauda_arvensis",
                  "Alces_alces","Anas_carolinensis",
                  "Anas_discors","Anas_platyrhynchos",
                  "Anthus_pratensis","Antilocapra_americana",
                  "Ardea_alba","Ateles_geoffroyi",
                  "Aythya_valisineria","Baeolophus_bicolor",
                  "Bupalus_piniaria","Buteo_lagopus",
                  "Callorhinus_ursinus","Canis_lupus",
                  "Capra_ibex","Ceratotherium_simum",
                  "Cervus_elaphus","Colaptes_auratus",
                  "Colinus_virginianus","Connochaetes_taurinus",
                  "Corvus_brachyrhynchos","Corvus_corone",
                  "Corvus_monedula","Cyanistes_caeruleus",
                  "Cygnus_olor","Dendrolimus_pini",
                  "Emberiza_schoeniclus","Engraulis_encrasicolus",
                  "Enhydra_lutris","Eriosoma_ulmi",
                  "Falco_peregrinus","Grus_americana",
                  "Gulosus_aristotelis","Haematopus_ostralegus",
                  "Istiophorus_platypterus","Larus_fuscus",
                  "Lepidorhombus_boscii","Lepidorhombus_whiffiagonis",
                  "Leptonychotes_weddellii","Lullula_arborea",
                  "Luscinia_megarhynchos","Lymantria_dispar",
                  "Macaca_sylvanus","Makaira_nigricans",
                  "Mallotus_villosus","Mareca_strepera",
                  "Marmota_sibirica","Melanerpes_erythrocephalus",
                  "Metopolophium_dirhodum","Myiarchus_crinitus",
                  "Myodes_gapperi","Napaeozapus_insignis",
                  "Neomonachus_schauinslandi","Odocoileus_virginianus",
                  "Oncorhynchus_gorbuscha","Operophtera_brumata",
                  "Orcinus_orca","Panthera_leo",
                  "Parus_major","Passerina_cyanea",
                  "Perdix_perdix","Phasianus_colchicus",
                  "Phoca_vitulina","Phylloscopus_collybita",
                  "Picoides_villosus","Procapra_gutturosa",
                  "Prunella_modularis","Sardina_pilchardus",
                  "Scomber_japonicus","Scomber_scombrus",
                  "Sterna_dougallii","Sturnus_vulgaris",
                  "Troglodytes_aedon","Urocitellus_parryii",
                  "Vanellus_vanellus","Vicugna_vicugna"))
obj_2 <- name.check(b,gpdd_p_core)

######################################################################

#### CONDUCT SIMPLE METRICS USING GRAPHS
#### Groups: Birds, Bony Fish, Insects, Mammals
library(ggplot2)

######################################################################
#### BOX 4
# get rid of repeating zeros
gpdd_fin$e[gpdd_fin$e=="#NAME?"]="-Inf"

# filter zeros
gpdd_zeros <- gpdd_fin %>%
  filter(as.numeric(e)<=1)

######################################################################
#### GRAPHS

gpdd_fin_yaya <- read.csv("final dataset.csv") # all
gpdd_p_yaya <- gpdd_fin_yaya %>%
  filter(p<=0.05) # trimmed nonpredictive dataset


#### DIMENSIONALITY - E

# get proportions of each class
gpdd_e <- gpdd_p_yaya %>%
  group_by(TaxonomicClass,E) %>%
  summarise(n=n()) %>%
  mutate(freq=n/sum(n))

# make E ordinal
gpdd_e$E <- ordered(gpdd_e$E)

# plot
windows()
gpdd_e %>%
  ggplot(aes(E,freq,fill=TaxonomicClass)) +
  geom_bar(stat="identity",position="dodge")+
  scale_x_discrete("Embedding Dimension (E)")+
  scale_y_continuous("Proportion")+
  scale_fill_manual(name="Taxonomic Class",
                    values=c("#56B4E9", "#009E73", "#F0E442", "#0072B2"))+
  theme_bw() +
  theme(#legend.position = "none",
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.title = element_text(size=20))

# do a boxplot as well by taxonomic class
fill <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2")

windows()
gpdd_p_yaya %>%
  ggplot(aes(x=TaxonomicClass,y=E))+
  geom_boxplot(fill=fill,position=position_dodge(width=0.82))+
  scale_x_discrete("Taxonomic Class")+
  scale_y_continuous("Embedding Dimension (E)")+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.title = element_blank())

# disaggregated by order
gpdd_z <- gpdd_p_yaya %>%
  group_by(TaxonomicOrder,E) %>%
  summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>% group_by(TaxonomicOrder) %>%
  mutate(total=sum(n))


windows()
gpdd_p_yaya %>%
  ggplot(aes(x=TaxonomicOrder,y=E))+
  geom_abline(intercept=0,slope=0,linetype="dashed",size=1)+
  geom_boxplot(position=position_dodge(width=0.82))+
  scale_x_discrete("Taxonomic Order")+
  scale_y_discrete("Embedding Dimension (E)")+
  theme_bw() +
  coord_flip() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.title = element_blank())

#### NONLINEARITY - THETA

gpdd_t <- gpdd_p_yaya %>%
  group_by(TaxonomicClass,nonlin) %>%
  summarise(n=n(),c=length(unique(MainID))) %>%
  mutate(freq=n/sum(n)) %>% group_by(TaxonomicClass) %>%
  mutate(total=sum(n))

# by class
windows()
gpdd_t %>%
  ggplot(aes(TaxonomicClass,
             freq,fill=nonlin)) +
  geom_bar(stat="identity",position="fill")+
  scale_x_discrete("Taxonomic Class")+
  scale_y_continuous("Proportion")+
  scale_fill_manual("",labels=c("Linear","Nonlinear"),
                    values=c("#0072B2", "#D55E00"))+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.title = element_text(size=20))

gpdd_z <- gpdd_p_yaya %>%
  group_by(TaxonomicOrder,nonlin) %>%
  summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>% group_by(TaxonomicOrder) %>%
  mutate(total=sum(n))

dumby_fxn <- function(x){
  ifelse(x$p>0.05,
         "a",
         ifelse(x$nonlin=="linear",
                "b",
                "c"))
}
gpdd_zz <- gpdd_fin_yaya %>%
  mutate(f=dumby_fxn(.data)) %>%
  group_by(TaxonomicOrder,f) %>%
  summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>% group_by(TaxonomicOrder) %>%
  mutate(total=sum(n))

gpdd_yy <- gpdd_fin_yaya %>%
  mutate(f=dumby_fxn(.data)) %>%
  group_by(TaxonomicClass,f) %>%
  summarise(n=n(),c=length(unique(MainID))) %>%
  mutate(freq=n/sum(n)) %>% group_by(TaxonomicClass) %>%
  mutate(total=sum(n))

# class with nonpredictable
windows()
gpdd_yy %>%
  #filter(total>=9) %>%
  # arrange(f) %>%
  ggplot(aes(TaxonomicClass,
             freq,fill=f)) +
  geom_bar(stat="identity",position="fill")+
  scale_x_discrete("Taxonomic Class")+
  scale_y_continuous("Proportion")+
  scale_fill_manual("",labels=c("Not Predictable","Linear","Nonlinear"),
                 values=c("#F0E442","#0072B2", "#D55E00"))+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.title = element_text(size=20))

# order with nonpredictable
windows()
gpdd_zz %>%
  filter(total>=9) %>%
  arrange(f) %>%
  ggplot(aes(TaxonomicOrder,
             freq,fill=f)) +
  geom_bar(stat="identity",position="fill")+
  scale_x_discrete("Taxonomic Order",
                     limits=rev(c("Hemiptera","Perciformes",
                                  "Pleuronectiformes","Rodentia",
                                  "Gadiformes","Anseriformes",
                                  "Charadriiformes","Falconiformes",
                                  "Lepidoptera","Passeriformes",
                                  "Artiodactyla","Carnivora")))+
  scale_y_continuous("Proportion")+
  scale_fill_manual("",labels=c("Not Predictable","Linear","Nonlinear"),
                    values=c("#F0E442","#0072B2", "#D55E00"))+
  coord_flip()+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.title = element_text(size=20))

# Forecast Ability - categorize by nonlin,lin, or not predictable
dumb_fxn <- function(x){
  ifelse(x$p>0.05,
  "Not Predictable",
  ifelse(x$nonlin=="linear",
  "Linear",
  "Nonlinear"))
}

gpdd_f <- gpdd_fin_yaya %>%
  mutate(f=dumb_fxn(.data))

# plot
windows()
gpdd_f %>%
  ggplot(aes(TaxonomicOrder,rho,fill=f))+
  geom_abline(intercept=0,slope=0,linetype="dashed",size=1)+
  geom_boxplot(position=position_dodge(width=0.82))+
  scale_x_discrete("Taxonomic Order")+
  scale_y_continuous("Forecast Skill (rho)")+
  theme_bw() +
  coord_flip() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.title = element_blank())

# dimensionality as a function of class + nonlinearity
windows()
gpdd_f %>%
  ggplot(aes(TaxonomicClass,E,fill=nonlin))+
  geom_abline(intercept=0,slope=0,linetype="dashed",size=1)+
  geom_boxplot(position=position_dodge(width=0.82))+
  scale_x_discrete("Taxonomic Class")+
  scale_y_continuous("Dimensionality (E)")+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.title = element_blank())

# Time Series Length as a function of Class + Nonlinearity
windows()
gpdd_p %>% ggplot(aes(TaxonomicClass,c,fill=nonlin))+
  geom_boxplot()

# time-series as examples of dimensionality and nonlinearity
# LOW E AND LINEAR
gpdd_woodcock <- gpdd_p %>% filter(MainID==10020)
gpdd_woodcock <- data.frame(gpdd_woodcock,seq(1,nrow(gpdd_woodcock)))
woodcock <- s_map(time_series=gpdd_woodcock$sc_Population,
                  lib=c(1,length(gpdd_woodcock$sc_Population)),
                  pred=c(1,length(gpdd_woodcock$sc_Population)),
                  E=4,theta=0.3,tp=1,stats_only = F)
woodcock_data <- woodcock$model_output[[1]]
gpdd_woodcock <- data.frame(gpdd_woodcock,woodcock_data)

windows()
gpdd_woodcock %>%
  ggplot(aes(x=seq.1..nrow.gpdd_woodcock..))+
  geom_abline(intercept=0,slope=0,linetype="dashed",size=1)+
  geom_line(aes(y=obs),color="gray",size=3)+
  geom_line(aes(y=pred),color="#0072B2",size=3)+
  scale_x_continuous("Time-steps")+
  scale_y_continuous("Standardized Abundance")+
  theme_bw() +
    theme(
      panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size=20),
      axis.text = element_text(size=16, color="black"),
      legend.text = element_text(size=16,color="black"),
      legend.title = element_blank())

# LOW E AND NONLINEAR
gpdd_vole <- gpdd %>% filter(MainID==10055) #10055
gpdd_vole <- data.frame(gpdd_vole,seq(1,nrow(gpdd_vole)))
vole <- s_map(time_series=gpdd_vole$sc_Population,
                  lib=c(1,length(gpdd_vole$sc_Population)),
                  pred=c(1,length(gpdd_vole$sc_Population)),
                  E=2,theta=1,tp=1,stats_only = F)
vole_data <- vole$model_output[[1]]
gpdd_vole <- data.frame(gpdd_vole,vole_data)

windows()
gpdd_vole %>%
  ggplot(aes(x=seq.1..nrow.gpdd_vole..))+
  geom_abline(intercept=0,slope=0,linetype="dashed",size=1)+
  geom_line(aes(y=obs),color="gray",size=3)+
  geom_line(aes(y=pred),color="#0072B2",size=3)+
  scale_x_continuous("Time-steps")+
  scale_y_continuous("Standardized Abundance")+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.title = element_blank())

# HIGH E AND LINEAR
gpdd_sole <- gpdd_p %>% filter(MainID==10063)
gpdd_sole <- data.frame(gpdd_sole,seq(1,nrow(gpdd_sole)))
sole <- s_map(time_series=gpdd_sole$sc_Population,
              lib=c(1,length(gpdd_sole$sc_Population)),
              pred=c(1,length(gpdd_sole$sc_Population)),
              E=7,theta=0,tp=1,stats_only = F)
sole_data <- sole$model_output[[1]]
gpdd_sole <- data.frame(gpdd_sole,sole_data)

windows()
gpdd_sole %>%
  ggplot(aes(x=seq.1..nrow.gpdd_sole..))+
  geom_abline(intercept=0,slope=0,linetype="dashed",size=1)+
  geom_line(aes(y=obs),color="gray",size=3)+
  geom_line(aes(y=pred),color="#0072B2",size=3)+
  scale_x_continuous("Time-steps")+
  scale_y_continuous("Standardized Abundance")+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.title = element_blank())

# HIGH E AND NONLINEAR
gpdd_aphid <- gpdd_p %>% filter(MainID==8648)
gpdd_aphid <- data.frame(gpdd_aphid,seq(1,nrow(gpdd_aphid)))
aphid <- s_map(time_series=gpdd_aphid$sc_Population,
              lib=c(1,length(gpdd_aphid$sc_Population)),
              pred=c(1,length(gpdd_aphid$sc_Population)),
              E=9,theta=2,tp=1,stats_only = F)
aphid_data <- aphid$model_output[[1]]
gpdd_aphid <- data.frame(gpdd_aphid,aphid_data)

windows()
gpdd_aphid %>%
  ggplot(aes(x=seq.1..nrow.gpdd_aphid..))+
  geom_abline(intercept=0,slope=0,linetype="dashed",size=1)+
  geom_line(aes(y=obs),color="gray",size=3)+
  geom_line(aes(y=pred),color="#0072B2",size=3)+
  scale_x_continuous("Time-steps")+
  scale_y_continuous("Standardized Abundance",c(-5,3))+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.title = element_blank())


####################################################################
#### ANALYSIS
#### GLM!!! - For Core Species
# Gdpp_fin = for rho; Gdpp_p = for E and nonlin

# lets scale data - just for fin
gpdd_fin_st <- gpdd_fin_core
gpdd_fin_st <- gpdd_fin_yaya

gpdd_fin_st <- gpdd_fin_st %>%
  mutate(f=dumb_fxn(.data))

gpdd_fin_st$c <- scale(gpdd_fin_st$c)
gpdd_fin_st$rho <- scale(gpdd_fin_st$rho)
gpdd_fin_st$E <- scale(gpdd_fin_st$E)
gpdd_fin_st$CV <- scale(gpdd_fin_st$CV)
gpdd_fin_st$time_mismatch <- scale(gpdd_fin_st$time_mismatch)
gpdd_fin_st <- data.frame(gpdd_fin_st,row.names = gpdd_fin_st$CorrectName)

# convinience wrapper
update_nested <- function(object, formula., ..., evaluate = TRUE){
  update(object = object, formula. = formula., data = object$model, ..., evaluate = evaluate)
}

#### Forecast Ability - RHO

# phylo try-out; remember to remove taxonomy
# remember that you have phyloPCA in
# and changed f to "gpdd_fin_stupid"

# try corPagel
rho_p_alt <- gls(rho ~ PC1+PC2 +
               c + f + CV + time_mismatch +
               SamplingUnits + SamplingUnits*CV,
             correlation = corPagel(1,a),
             data=gpdd_fin_st, method="ML")

rho_p <- gls(rho ~ PC1.1+PC2.1 +
               c + f + CV + time_mismatch +
               SamplingUnits + SamplingUnits*CV,
             correlation = corBrownian(phy=a),
             data=gpdd_fin_st, method="ML")


rho_1 <- glm(rho ~ TaxonomicClass + PC1.1 + PC2.1 +
               c + f + CV + time_mismatch +
               SamplingUnits + SamplingUnits*CV,
               data=gpdd_fin_st)

# model selection
rho_drop1 <- update_nested(rho_1, .~. -CV:SamplingUnits)
rho_drop2 <- update_nested(rho_drop1,.~.-SamplingUnits)
rho_drop3 <- update_nested(rho_drop2,.~.-c)
rho_drop4 <- update_nested(rho_drop3,.~.-PC2)
rho_drop5 <- update_nested(rho_drop4,.~.-correct.value)
rho_drop6 <- update_nested(rho_drop5,.~.-time_mismatch)
rho_drop7 <- update_nested(rho_drop6,.~.-TaxonomicClass)

#### Dimensionality - E

gpdd_p_st <- gpdd_p_core
gpdd_p_st <- gpdd_p_yaya
gpdd_p_st <- gpdd_p

# random stuff for reviewer
#gpdd_p_st$minageN <- gpdd_p_st$MinAge/gpdd_p_st$c
#gpdd_p_st$longevityN <- gpdd_p_st$Lifesp/gpdd_p_st$c

gpdd_p_st$TrL <- factor(gpdd_p_st$TrL)
gpdd_p_st$c <- scale(gpdd_p_st$c)
gpdd_p_st$rho <- scale(gpdd_p_st$rho)
gpdd_p_st$CV <- scale(gpdd_p_st$CV)
gpdd_p_st$time_mismatch <- scale(gpdd_p_st$time_mismatch)
gpdd_p_st$minageN <- scale(gpdd_p_st$minageN)
gpdd_p_st$longevityN <- scale(gpdd_p_st$longevityN)

#library(ordinal) - using ordinal is weird/messy??
library(ordinal)

gpdd_p_st$E <- ordered(gpdd_p_st$E)

dim_p_alt <- gls(E~PC1.1 + PC2.1 + TrL  +
               time_mismatch + c + CV + SamplingUnits,
             correlation = corPagel(seq(0,1,0.1),phy=b,fixed=T),
             data=gpdd_p_st, method="ML")

dim_p <- gls(E~PC1.1 + PC2.1 + TrL  +
               time_mismatch + c + CV + SamplingUnits,
             correlation = corBrownian(phy=b),
             data=gpdd_p_st, method="ML")

dim_1 <- clm(E ~ TaxonomicClass + PC1 + PC2 + TrL  +
             time_mismatch + c + CV + SamplingUnits,
             data=gpdd_p_st)

# stuff for reviewers


dim_whatev <- clm(E ~ longevityN,data=gpdd_p_st)

# summary
dim_drop1 <- update_nested(dim_1, .~. -SamplingUnits)
dim_drop2 <- update_nested(dim_drop1, .~. -TrL)
dim_drop3 <- update_nested(dim_drop2, .~. -TaxonomicClass)
dim_drop4 <- update_nested(dim_drop3, .~. -time_mismatch)
dim_drop5 <- update_nested(dim_drop4, .~. -CV)
dim_drop6 <- update_nested(dim_drop5, .~. -PC2)

#### Nonlinearity - Theta
library(arm)
library(phylolm)

# convert to numerals
gpdd_p_st$nonlin <- as.character(gpdd_p_st$nonlin) #extra
gpdd_p_st$nonlin[gpdd_p_st$nonlin=="linear"]=0
gpdd_p_st$nonlin[gpdd_p_st$nonlin=="nonlinear"]=1
gpdd_p_st$nonlin <- as.numeric(gpdd_p_st$nonlin)

ifelse(gpdd_p_st$nonlin[gpdd_p_st$nonlin=="linear"],
       gpdd_p_st$nonlin[gpdd_p_st$nonlin] <- 0,
       gpdd_p_st$nonlin[gpdd_p_st$nonlin] <- 1)

# insect stuff for reviewer
insects <- gpdd_p_st %>% filter(TaxonomicClass=="Insecta")


theta_p <- phyloglm(nonlin~ PC1.1 + PC2.1 +
                      E + rho + c + CV +
                      time_mismatch ,
                    method="logistic_IG10",phy=b,
                    data=gpdd_p_st)

theta_p_alt <- binaryPGLMM(nonlin~ PC1.1 + PC2.1 +
                             E + rho + c + CV +
                             time_mismatch,phy=b,
                           data=gpdd_p_st)

theta_1 <- glm(nonlin ~ PC1.1 + PC2.1 +
                 E + rho + c + CV+
                 time_mismatch,
               data=gpdd_p_st, family="binomial")


# stuff for reviewers
#blahstuff <- glm(nonlin~longevityN, data=gpdd_p_st,family="binomial")

theta_1 <- bayesglm(nonlin ~ TaxonomicClass + PC1 + PC2 +
               E + rho + c + CV + SamplingUnits +
               time_mismatch + SamplingUnits*CV,
               data=gpdd_p_st, family="binomial")

insects_1 <- bayesglm(nonlin ~ PC1 + PC2 + E + rho + c + CV+
                        SamplingUnits + time_mismatch +
                        SamplingUnits*CV,
                      data=insects, family="binomial")

insects_2 <- update_nested(insects_1, .~. -SamplingUnits:CV)
insects_3 <- update_nested(insects_2, .~. -SamplingUnits)
insects_4 <- update_nested(insects_3, .~. -c)
insects_5 <- update_nested(insects_4, .~. -PC1)
insects_6 <- update_nested(insects_5,.~. -CV)
insects_7 <- update_nested(insects_6,.~. -PC2)
insects_8 <- update_nested(insects_7,.~. -time_mismatch)
insects_9 <- update_nested(insects_8,.~. -E)

# model selection
theta_drop1 <- update_nested(theta_1, .~. -SamplingUnits:CV)
theta_drop2 <- update_nested(theta_drop1, .~. -SamplingUnits)
theta_drop3 <- update_nested(theta_drop2, .~. -time_mismatch)
theta_drop4 <- update_nested(theta_drop3, .~. -PC1)
theta_drop5 <- update_nested(theta_drop4, .~. -PC2)
theta_drop6 <- update_nested(theta_drop5, .~. -CV)

# graph of the interaction
# why are long-lived insects so nonlinear???

windows()
plot(seq(-1,4,length.out=10),seq(-50,50,length.out = 10),type="n")
lines(seq(-1,4,length.out=10),-8.549e-01+(8.644e-01*seq(-1,8,1)),
      col="red")#birds
lines(seq(-1,4,length.out=10),7.9851+(11.8544*seq(-1,8,1)),
      col="green")#ins
lines(seq(-1,4,length.out=10),0.3631+(-0.1011*seq(-1,8,1)),
      col="purple")#mam
lines(seq(-1,4,length.out=10),-0.0415+(0.1873*seq(-1,8,1)),
      col="blue")#fish


# with mammals
gpdd_mam <- gpdd_p %>% filter(TaxonomicClass=="Mammalia")

theta_mam <- glm(nonlin ~ Len + MinAge +
                 Lifesp + Fert + c + SamplingUnits +
                   SamplingUnits*cent_sd, data=gpdd_mam, family="binomial")



#####################################################################
#### OTHER RANDOM STUFF

write.csv(gpdd_4 %>% filter(TaxonomicClass=="Mammalia") %>% group_by(TaxonName) %>% summarize(n()),"mammal_list.csv")
write.csv(gpdd_4 %>% filter(TaxonomicClass=="Aves") %>% group_by(TaxonName) %>% summarize(n()),"bird_list.csv")

gpdd_filter %>% filter(TaxonomicClass == "Osteichthyes") %>%
  group_by(TaxonName) %>%
     summarize(Length=length(unique(MainID)),
               Species=length(unique(TaxonName)))

# create csv of all data used - remember! use this for graphs!
gpdd_fin_yaya <- gpdd_fin %>% group_by(MainID) %>%
  sample_n(1)
gpdd_p_yaya <- gpdd_p %>% group_by(MainID) %>%
  sample_n(1)
write.csv(gpdd_fin_yaya,"final dataset.csv")

  summarize_(SamplingUnits=SamplingUnits,Notes.x=Notes.x,
            TaxonName=TaxonName,CommonName=CommonName,
            TaxonomicClass=TaxonomicClass,TaxonomicOrder=TaxonomicOrder,
            TaxonomicFamily=TaxonomicFamily,TaxonomicGenus=TaxonomicGenus,
            Mass=Mass,Len=Len,Wing=Wing,MinAge=MinAge,Lifesp=Lifesp,
            Fert=Fert,TrL=TrL,N=c,CV=CV,E=E,G=time_mismatch,
            nonlin=nonlin,rho=rho,predictable=f)

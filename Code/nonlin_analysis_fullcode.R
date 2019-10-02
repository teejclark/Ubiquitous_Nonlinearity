#### Ubiquitous Nonlinearity Code
#### T.J. Clark and Angela Luis - 10/2/19

rm(list=ls(all=T))

library("rgpdd")
library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")

###################################################################

#### 1) DATA PREPARATION
# most of the procedure to bring together/create the dataset
# can skip by just importing full dataset stored in Github

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

# Join with other data
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

# Attach Corrected Data - IMPORTANT
# these are grey datasets that didn't match the GPDD format
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
#### 2) DATA ANALYSIS
# create functions and run EDM on dataset

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
#### Box A
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

# or load random phase matrix created using
# phase-randomized surrogate series, Ebisuzaki (1997)

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
#### 3) SUMMARY STATS
# looking at the data

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


# PCA ANALYSIS 
# to bring different life history traits into 2 variables
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
# in effect, sampling.
gpdd_fin_core <- gpdd_fin %>%
  group_by(TaxonName) %>%
  sample_n(1)


####################################################################
#### 4) ANALYSIS

#### GLM!!! - For Core Species
# Gdpp_fin = for rho; Gdpp_p = for E and nonlin

# lets scale data - just for fin
gpdd_fin_st <- gpdd_fin_core
gpdd_fin_st <- gpdd_fin_yaya

gpdd_fin_st <- gpdd_fin_st %>%
  mutate(f=dumb_fxn(.data))

# scale predictors
gpdd_fin_st$c <- scale(gpdd_fin_st$c)
gpdd_fin_st$rho <- scale(gpdd_fin_st$rho)
gpdd_fin_st$E <- scale(gpdd_fin_st$E)
gpdd_fin_st$CV <- scale(gpdd_fin_st$CV)
gpdd_fin_st$time_mismatch <- scale(gpdd_fin_st$time_mismatch)
gpdd_fin_st <- data.frame(gpdd_fin_st,row.names = gpdd_fin_st$CorrectName)

# convinience wrapper for model selection
update_nested <- function(object, formula., ..., evaluate = TRUE){
  update(object = object, formula. = formula., data = object$model, ..., evaluate = evaluate)
}

#### Forecast Ability - RHO

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

gpdd_p_st$TrL <- factor(gpdd_p_st$TrL)
gpdd_p_st$c <- scale(gpdd_p_st$c)
gpdd_p_st$rho <- scale(gpdd_p_st$rho)
gpdd_p_st$CV <- scale(gpdd_p_st$CV)
gpdd_p_st$time_mismatch <- scale(gpdd_p_st$time_mismatch)
gpdd_p_st$minageN <- scale(gpdd_p_st$minageN)
gpdd_p_st$longevityN <- scale(gpdd_p_st$longevityN)

library(ordinal)

gpdd_p_st$E <- ordered(gpdd_p_st$E)

dim_1 <- clm(E ~ TaxonomicClass + PC1 + PC2 + TrL  +
             time_mismatch + c + CV + SamplingUnits,
             data=gpdd_p_st)

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

theta_1 <- bayesglm(nonlin ~ TaxonomicClass + PC1 + PC2 +
               E + rho + c + CV + SamplingUnits +
               time_mismatch + SamplingUnits*CV,
               data=gpdd_p_st, family="binomial")

# model selection
theta_drop1 <- update_nested(theta_1, .~. -SamplingUnits:CV)
theta_drop2 <- update_nested(theta_drop1, .~. -SamplingUnits)
theta_drop3 <- update_nested(theta_drop2, .~. -time_mismatch)
theta_drop4 <- update_nested(theta_drop3, .~. -PC1)
theta_drop5 <- update_nested(theta_drop4, .~. -PC2)
theta_drop6 <- update_nested(theta_drop5, .~. -CV)


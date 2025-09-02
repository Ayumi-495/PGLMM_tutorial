######## MCMCglmm - for non-MCA analyses (VS only, VS + AC, Diurnal only) ########



library(MCMCglmm)
library(phangorn)
library(ape)
source("C:/Users/ve21130/OneDrive - University of Bristol/Thesis/Post-thesis/New manuscript/To submit/Scripts/VIF_functions.R")   # Include custom functions
library(caper)
library(stableGR)
library(mcmcse)
library(brio)


# read in data - change file location as necessary
# If looking at male colouration, use "Primate_data_male.csv", for female colouration, "Primate_data_female.csv"
primateData <- read.csv("C:/Users/ve21130/OneDrive - University of Bristol/Thesis/Post-thesis/New manuscript/To submit/Datasets/Primate_data_male.csv")
data <- primateData  


###### crop data ######

# Remove rows with missing data
data <- subset(data, !is.na(vs_male))     # omit NAs by VS predictor
data <- subset(data, !is.na(activity_cycle)) # omit NAs by AC
data <- subset(data, !is.na(redpeachpink_facial_skin))  # omit NAs by response variable
data <- subset(data, !is.na(social_group_size))  # omit NAs by response variable
data <- subset(data, !is.na(multilevel))  # omit NAs by response variable

# Rename activity cycle levels to make diurnal the baseline (for all analyses except diurnal-only)
data$activity_cycle[data$activity_cycle == "di"] <- "a_di"






#a random 100 topologies of the Upham et al. 2019 mammal tree

t100<-read.nexus("C:/Users/ve21130/OneDrive - University of Bristol/Thesis/Post-thesis/New manuscript/To submit/Scripts/trees100m.nex") # Random 100 of 10,000 topologies

#trim to the tree to match the dataset
tree<-t100[[1]]

t100<-lapply(t100,drop.tip,tip=setdiff(tree$tip.label,data$PhyloName)) #trim out everything from the tree that's not in the dataset

tree<-t100[[1]] #select one tree for trimming purposes



################################################################################################################
############################ Presence/absence of colour (logistic regression)  #####################
################################################################################################################


######
#set up a dummy run
######
i=1 #this is arbitrary

tree<-t100[[i]]  # pull a tree from the distribution (number i - starts at one)
tree<-nnls.tree(cophenetic(tree),tree,method = "ultrametric") # force tree to be ultrametric - all tips equidistant from root

animalA<-inverseA(tree)$Ainv # invert covariance matrix for use by MCMCglmm 

#prior for a logistic regression [discrete dependent variables]
# [Gives the prior covariance matrix for fixed effects]
gelmanprior<-list(B=list(mu=rep(0,9), #set mu up for the number of coefficients you're estimating -- one for intercept, one for slope of each predictor variable (N.B. for categorical variables, each level counts as a predictor)
                         # For VS+AC analyses - hash out if running VS+AC models                         
                         V=gelman.prior(~ vs_male + activity_cycle + multilevel + social_group_size + vs_male*social_group_size, #estimate a prior based on the equation that you're running
                                        data = data,  scale=1+pi^2/3)), 
                  R = list(V = 1,fix = 1), G=list(G1 = list(V = 1E-10, nu = -1 )))    # Based on Hadfield


######## Dummy run ###########

# For VS+AC analyses - hash out for VS + AC analyses
mod<-MCMCglmm(redpeachpink_facial_skin ~ vs_male + activity_cycle + multilevel + social_group_size + vs_male*social_group_size, 
              random=~PhyloName,                  # [Sets phylogeny as a random factor]
              ginverse=list(PhyloName=animalA), # [list of sparse inverse matrices proportional to covariance structure of the random effects - name nees to correspond to random effects column in data]
              prior = gelmanprior,   
              verbose=TRUE,       # [Prints MH diagnostics to screen]
              family="categorical", #family for logistic regression -- use "gaussian" if you're running a linear regression
              data = data,
              nitt=85000,   # [number of MCMC iterations]
              thin=80,      # [thinning interval - reduces memory use, and helps to reduce autocorrelation of sample and increase independence]
              burnin=5000,  # [discard unconverged values]
              pl=TRUE,      # [save posterior distribution of latent (inferred) variables]
              pr=TRUE,      # [save posterior distribution of random variables]
              slice=TRUE) 

Final.mod<-mod #set up a structure that we'll populate with the real model
Final.mod$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,] # [VCV is posterior distrib of covariance matrices]
Final.mod$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] # [Sol is posterior distrib of MME solutions - includes fixed effects]
Final.mod$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] # [Liab is posterior distrib of latent variables]


nsamp.l<-nrow(mod$VCV)   
start1.l=list(R=mod$VCV[nsamp.l,"units"], G=list(G1=mod$VCV[nsamp.l,"PhyloName"])) 

# Create model file to save
save(Final.mod,file="80k-col-angol-VS+AC+ML+SGS+int-M-redpeachpink-facial-skin-VS-MA.Rdata")





############### Loop over 100 trees ###################


for(i in 1:100){ #loop through 100 trees
  tree<-t100[[i]] #select the ith tree 
  tree<-nnls.tree(cophenetic(tree),tree, method = "ultrametric") #the tree is non-ultrametric, which is annoying, so we force it to be ultrametric
  
  
  animalA<-inverseA(tree)$Ainv 
  
  # For VS+AC analyses - hash out for VS + AC analyses
  mod<-MCMCglmm(redpeachpink_facial_skin ~ vs_male + activity_cycle + multilevel + social_group_size + vs_male*social_group_size,  
                random = ~PhyloName,
                ginverse=list(PhyloName=animalA), 
                prior = gelmanprior, 
                verbose=FALSE, 
                family="categorical", 
                start= start1.l,
                data = data,
                nitt=80000,   # [number of MCMC iterations]
                thin=7500,      # [thinning interval - reduces memory use, and maybe helps to reduce autocorrelation of sample and increase independence]
                burnin=5000,
                pl=TRUE,
                pr=TRUE,
                slice=TRUE)
  
  print(i) #print which tree you're on (for your sanity as the loop runs)
  
  Final.mod$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,] # [VCV is posterior distrib of covariance matrices]
  Final.mod$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] # [Sol is posterior distrib of MME solutions (???) - includes fixed effects]
  Final.mod$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] # [Liab is posterior distrib of latent variables]
  
  nsamp.l<-nrow(mod$VCV)
  start1.l=list(R=mod$VCV[nsamp.l,"units"], G=list(G1=mod$VCV[nsamp.l,"PhyloName"]))
  
  if(i > 99){
      save(Final.mod,file="80k-col-angol-VS+AC+ML+SGS+int-M-redpeachpink-facial-skin-VS-MA.Rdata") #save model file at end of run
  }
}



# Results

# Plot trace/density plots
plot(Final.mod)

# See summary of results
summary(Final.mod)

# Compile results into table and copy to clipboard
sum1 <- as.data.frame(summary(Final.mod)$solutions)
write.excel(sum1, row.names = TRUE)

# Copy DIC to clipboard
write.excel(summary(Final.mod)$DIC)

# See standard deviations
summary(Final.mod$Sol[,1:3])



######## Diagnostics  #########

### Autocorrelation

autocorrMCMCglmm(Final.mod, 9)


autocorr(Final.mod$Sol[,1:9])  # Check for convergence in the fixed effects
autocorr(Final.mod$VCV)    # Check for convergence in the random effects (0.1 is a good threshold)


# Plot autocorrelation
autocorr.plot(Final.mod$VCV)    # Check for convergence in the random effects (0.1 is a good threshold)
autocorr.plot(Final.mod$Sol[,1:9])  # Check for convergence in the fixed effects

# Use Heidelberger & Welch test  (stationary/halfwidth)
heidel.diag(Final.mod$VCV)
heidel.diag(Final.mod$Sol[, 1:9])

# Geweke diagnostic (want value between -2 and 2)
geweke.diag(Final.mod$VCV)
geweke.diag(Final.mod$Sol[, 1:9])
geweke.plot(Final.mod$Sol[, 1:9])

### Variance inflation factors (to check problems with multicollinearity)

vif.MCMCglmm(Final.mod)    # Check VIFs - Austin Frank method - N.B. this isn't appropriate for categorical predictors - it interprets it as multiple dummy variables

# This computes GVIFS - appropriate for categorical predictor variables
VIFdata <- cbind(data$vs_male, data$activity_cycle, data$multilevel, as.numeric(data$social_group_size))  # Create dataset of predictor variables only to use with corvif function - change predictor variables as necessary
colnames(VIFdata) <- c("vs_male", "activity_cycle", "multilevel", "social_group_size")
corvif(VIFdata)   # evaluates GVIFs



# For checking convergence by comparing multiple runs, if required - computes Gelman diagnostic

run1 <- Final.mod
run2 <- Final.mod

geldiag <- gelman.diag(mcmc.list(run1$Sol[,1:9], run2$Sol[,1:9]))   # Check for convergence

stable.GR(Final.mod$Sol[,1:9])   #  Calculates improved GR diagnostic (Vats & Knudson 2021) - only need 1 chain

n.eff(Final.mod$Sol[,1:9])


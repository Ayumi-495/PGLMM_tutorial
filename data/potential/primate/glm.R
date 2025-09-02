# Non-phylogenetic tests using glm() function

source("C:/Users/ve21130/OneDrive - University of Bristol/Thesis/Post-thesis/New manuscript/To submit/Scripts/VIF_functions.R")   # Include custom functions

# read in data - change file location as necessary
primateData <- read.csv("C:/Users/ve21130/OneDrive - University of Bristol/Thesis/Post-thesis/New manuscript/To submit/Datasets/Primate_data_male.csv")



################################################################################
############### Including AC, SGS & ML as covariate - no interaction term ######
################################################################################

# Remove rows with missing data

data <- primateData 
 
data <- subset(data, !is.na(vs_male))     # omit NAs by VS predictor
data <- subset(data, !is.na(red_pelage_body_limbs))  # omit NAs by response variable
data <- subset(data, !is.na(activity_cycle)) # omit NAs by AC
data <- subset(data, !is.na(social_group_size)) # omit NAs by social group size
data <- subset(data, !is.na(multilevel)) # omit NAs by multilevel y/n

# Rename activity cycle levels to make diurnal the baseline 
data$activity_cycle[data$activity_cycle == "di"] <- "a_di"


# Run regression

# Basic (w/ no confounders)
basemod_h <- glm(red_pelage_head ~ vs_male + activity_cycle + multilevel + social_group_size, family = binomial(), data = data)
basemod_bl <- glm(red_pelage_body_limbs ~ vs_male + activity_cycle + multilevel + social_group_size, family = binomial(), data = data)
basemod_t <- glm(red_pelage_tail ~ vs_male + activity_cycle + multilevel + social_group_size, family = binomial(), data = data)



# Compile results into table and copy to clipboard
sum1 <- as.data.frame(summary(basemod_h)$coefficients)
write.excel(sum1, row.names = TRUE, col.names = FALSE)

sum1 <- as.data.frame(summary(basemod_bl)$coefficients)
write.excel(sum1, row.names = TRUE, col.names = FALSE)

sum1 <- as.data.frame(summary(basemod_t)$coefficients)
write.excel(sum1, row.names = TRUE, col.names = FALSE)


# Same, but for sexual skin colouration
data <- primateData  

data <- subset(data, !is.na(vs_male))     # omit NAs by VS predictor
data <- subset(data, !is.na(comp_red_gen_female))  # omit NAs by response variable
data <- subset(data, !is.na(activity_cycle)) # omit NAs by AC
data <- subset(data, !is.na(social_group_size)) # omit NAs by social group size
data <- subset(data, !is.na(multilevel)) # omit NAs by multilevel y/n

# Rename activity cycle levels to make diurnal the baseline 
data$activity_cycle[data$activity_cycle == "di"] <- "a_di"

# 
basemod_ss <- glm(comp_red_gen_female ~ vs_male + activity_cycle + multilevel + social_group_size, family = binomial(), data = data)


# Compile results into table and copy to clipboard
sum1 <- as.data.frame(summary(basemod_ss)$coefficients)
write.excel(sum1, row.names = TRUE, col.names = FALSE)


# And for facial skin colouration
data <- primateData  

data <- subset(data, !is.na(vs_male))     # omit NAs by VS predictor
data <- subset(data, !is.na(redpeachpink_facial_skin))  # omit NAs by response variable
data <- subset(data, !is.na(activity_cycle)) # omit NAs by AC
data <- subset(data, !is.na(social_group_size)) # omit NAs by social group size
data <- subset(data, !is.na(multilevel)) # omit NAs by multilevel y/n

# Rename activity cycle levels to make diurnal the baseline 
data$activity_cycle[data$activity_cycle == "di"] <- "a_di"

# Basic (w/ no confounders)
basemod_rppfs <- glm(redpeachpink_facial_skin ~ vs_male + activity_cycle + multilevel + social_group_size, family = binomial(), data = data)


# Compile results into table and copy to clipboard
sum1 <- as.data.frame(summary(basemod_rppfs)$coefficients)
write.excel(sum1, row.names = TRUE, col.names = FALSE)

# Get aic scores for all models
summary(basemod_h)$aic
summary(basemod_bl)$aic
summary(basemod_t)$aic
summary(basemod_ss)$aic
summary(basemod_rppfs)$aic






################################################################################
############### Including AC, SGS & ML as covariate - with interaction #########
################################################################################


# Remove rows with missing data
data <- primateData 
 
data <- subset(data, !is.na(vs_male))     # omit NAs by VS predictor
data <- subset(data, !is.na(red_pelage_body_limbs))  # omit NAs by response variable
data <- subset(data, !is.na(activity_cycle)) # omit NAs by AC
data <- subset(data, !is.na(social_group_size)) # omit NAs by social group size
data <- subset(data, !is.na(multilevel)) # omit NAs by multilevel y/n

# Rename activity cycle levels to make diurnal the baseline 
data$activity_cycle[data$activity_cycle == "di"] <- "a_di"

# Check VIFs
VIFdata <- cbind(data$vs_male, data$activity_cycle, as.numeric(data$social_group_size), as.factor(data$multilevel))  # Create dataset of predictor variables only to use with corvif function - change predictor variables as necessary
corvif(VIFdata)

# Run regression

# Basic (w/ no confounders)
basemod_h <- glm(red_pelage_head ~ vs_male + activity_cycle + multilevel + social_group_size + vs_male*social_group_size, family = binomial(), data = data)
basemod_bl <- glm(red_pelage_body_limbs ~ vs_male + activity_cycle + multilevel + social_group_size + vs_male*social_group_size, family = binomial(), data = data)
basemod_t <- glm(red_pelage_tail ~ vs_male + activity_cycle + multilevel + social_group_size + vs_male*social_group_size, family = binomial(), data = data)


# Compile results into table and copy to clipboard
sum1 <- as.data.frame(summary(basemod_h)$coefficients)
write.excel(sum1, row.names = TRUE, col.names = FALSE)

sum1 <- as.data.frame(summary(basemod_bl)$coefficients)
write.excel(sum1, row.names = TRUE, col.names = FALSE)

sum1 <- as.data.frame(summary(basemod_t)$coefficients)
write.excel(sum1, row.names = TRUE, col.names = FALSE)


# Same, but for sexual skin colouration
data <- primateData 

 
data <- subset(data, !is.na(vs_male))     # omit NAs by VS predictor
data <- subset(data, !is.na(comp_red_gen_female))  # omit NAs by response variable
data <- subset(data, !is.na(activity_cycle)) # omit NAs by AC
data <- subset(data, !is.na(social_group_size)) # omit NAs by social group size
data <- subset(data, !is.na(multilevel)) # omit NAs by multilevel y/n

# Rename activity cycle levels to make diurnal the baseline 
data$activity_cycle[data$activity_cycle == "di"] <- "a_di"

# 
basemod_ss <- glm(comp_red_gen_female ~ vs_male + activity_cycle + multilevel + social_group_size + vs_male*social_group_size, family = binomial(), data = data)


# Compile results into table and copy to clipboard
sum1 <- as.data.frame(summary(basemod_ss)$coefficients)
write.excel(sum1, row.names = TRUE, col.names = FALSE)


# And for facial skin colouration
data <- primateData 

 
data <- subset(data, !is.na(vs_male))     # omit NAs by VS predictor
data <- subset(data, !is.na(redpeachpink_facial_skin))  # omit NAs by response variable
data <- subset(data, !is.na(activity_cycle)) # omit NAs by AC
data <- subset(data, !is.na(social_group_size)) # omit NAs by social group size
data <- subset(data, !is.na(multilevel)) # omit NAs by multilevel y/n

# Rename activity cycle levels to make diurnal the baseline 
data$activity_cycle[data$activity_cycle == "di"] <- "a_di"

# Basic (w/ no confounders)
basemod_rppfs <- glm(redpeachpink_facial_skin ~ vs_male + activity_cycle + multilevel + social_group_size + vs_male*social_group_size, family = binomial(), data = data)


# Compile results into table and copy to clipboard
sum1 <- as.data.frame(summary(basemod_rppfs)$coefficients)
write.excel(sum1, row.names = TRUE, col.names = FALSE)

# Get aic scores for all models
summary(basemod_h)$aic
summary(basemod_bl)$aic
summary(basemod_t)$aic
summary(basemod_ss)$aic
summary(basemod_rppfs)$aic





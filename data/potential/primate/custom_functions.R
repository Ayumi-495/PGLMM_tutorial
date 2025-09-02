# VIF & lambda estimation functions - don't mess with these


############### Function allowing estimation of lambda ####################
# Based on calculation from de Villemereuil and Nakagawa 2014 (from Garamszegi 2014 Modern phylogenetic comparative methods - Ch 11.2.1)
# Essentially this computes how much of the variance is explained by the random effects, which in this case is the phylogeny

lamCalc <- function(mod, phylo){
  phylo <- as.character(phylo)
  lambda <- mod$VCV[, phylo] / (mod$VCV[, phylo] + mod$VCV[, "units"])
  return(lambda)
}

#################### VIF function ########################### (Austin Frank)

vif.MCMCglmm <- function (fit, intercept.columns = c(1)) {
  nF <- fit$Fixed$nfl
  v <- cov(as.matrix(fit$X[,1:nF]))
  nam <- colnames(fit$Sol[,1:nF])
  v <- v[-intercept.columns, -intercept.columns, drop = FALSE]
  nam <- nam[-intercept.columns]
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}


########### VIF function for glm models  (Robert MacDonald) ###########################
########### N.B. for this to work, when call glm() function need to specify x = TRUE

glmVIF <- function (fit, intercept.columns = c(1)) {
  if ("x" %in% names(fit) ==  FALSE) stop("The model does not contain a model matrix. Suggest rerunning model with x = TRUE")
  nF <- fit$rank
  v <- cov(as.matrix(fit$x[,1:nF]))
  nam <- rownames(data.frame(testmod$effects[1:nF]))
  v <- v[-intercept.columns, -intercept.columns, drop = FALSE]
  nam <- nam[-intercept.columns]
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

############### Alternative VIF function ###################

#Library files for courses provided by: Highland Statistics Ltd.
#To cite these functions, use:
#Mixed effects models and extensions in ecology with R. (2009).
#Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.

#Copyright Highland Statistics LTD.

#####################################################################
#VIF FUNCTION.
#To use:  corvif(YourDataFile)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}


#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}
#END VIF FUNCTIONS



# Little function to copy to clipboard so can be pasted to Excel
write.excel <- function(x,row.names=FALSE,col.names=FALSE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}


# Function to check autocorrelation in first non-zero lag of MCMCglmm fixed effects
autocorrMCMCglmm <- function(mod, var) {
  var2 <- var*var
  vals <- data.frame(coda::autocorr(mod$Sol[,1:var]))
  count <- 0
  for(i in 1:var2){
    if(vals[2, i] >= 0.1){
      return(paste0("Uh oh, problem in column <", colnames(vals[i]), ">; Lag = ", vals[2, i]))
      count <- count + 1
    }
  }
  if(count == 0){
    return("All good!")
  }
}

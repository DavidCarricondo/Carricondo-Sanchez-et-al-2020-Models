
#JAGS specification of the models testing the behavior of wolves in relation to 
#human infrastructures in central and southeastern Scandinavia from 2001 to 2017.


#Conditional logistic model

model {
  for (i in 1:nStrata) {
    for (j in 1:nSteps) {
      
      log(phi[i,j]) <- beta_distb * cdist_build[i,j] + 
        beta_dens * cdensbuild[i,j] +
        beta_distm * cdist_mainroad[i,j]+
        beta_dists * cdist_secroad[i,j]+
        beta_for * for_edgekm[i,j]
      
      p[i,j] <- phi[i,j]/sum(phi[i,1:6]) 
    }
    
    y[i,1:nSteps] ~ dmulti(p[i,1:nSteps],1)       
  }
  
  beta_distb ~ dnorm(0,5)
  beta_dens ~ dnorm(0,5)
  beta_distm ~ dnorm(0,5)
  beta_dists ~ dnorm(0,5)
  beta_for ~ dnorm(0,5)
}

#Weighted mixed effects model with covariates:
  
  model{
    #Priors
    
    tau <- pow(sd, -2)
    sd ~ dunif(0, 10)
    mean.alpha ~ dnorm(0, 5)
    tau.alpha <- pow(sd.alpha, -2)
    sd.alpha ~ dunif(0, 10) # sd hyperparam
    
    for (i in 1:n.par){
      beta[i] ~ dnorm(0,5)}
    
    for (j in 1:n.ind){
      alpha[j] ~ dnorm(mean.alpha, tau.alpha)
    }
    
    #likelihood
    
    for (i in 1:n.obs){
      y[i] ~ dnorm(mean[i], tau*wt[i]) ##Wt is already inversed in the data specification wt=1/se
      mean[i] <- alpha[id[i]] + inprod(x[i, 1:n.par],beta[1:n.par])
    }
    #Derived quantities
    
    for (i in 1:n.obs){
      resi[i] <- y[i]-mean[i]  
      log_lik[i] <- logdensity.norm(y[i], mean[i], tau*wt[i])  #log likelihood calculated with the model weights
    }
    
  }


#Weighted mixed effects model without covariates (null models):
  
  model{
    #Priors
    tau <- pow(sd, -2)
    sd ~ dunif(0, 10)
    mean.alpha ~ dnorm(0, 5)
    tau.alpha <- pow(sd.alpha, -2)
    sd.alpha ~ dunif(0, 10) # sd hyperparam
    
    for (j in 1:n.ind){
      alpha[j] ~ dnorm(mean.alpha, tau.alpha)
    }
    
    #likelihood
    
    for (i in 1:n.obs){
      y[i] ~ dnorm(mean[i], tau*wt[i])  ##Wt is already inversed in the data specification wt=1/se
      mean[i] <- alpha[id[i]]
    }
    
    #Derived quantities
    
    for (i in 1:n.obs){
      resi[i] <- y[i]-mean[i]  
      log_lik[i] <- logdensity.norm(y[i], mean[i], tau*wt[i])  #log likelihood calculated with the model weights
    }
  }

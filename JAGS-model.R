model{  
  ####################################
  # Indices and parameter descriptions
  ####################################
  # lmu.N log-transformed abundance index after effort (effort.effect), circle, and distance corrections (dist offset)
  # r= population growth rate on log scale, 
  # lambda= backtransformed population growth rate
  # sig.obs= observation error of counts, sig.proc= process error of abundanceindices
  # k= count in ncounts, y= year (also yr), 
  # s= state (also str) in nstrata ,
  
  # Priors
  mean.r ~ dnorm(0, 0.001)
  mean.lambda<- exp(mean.r)
  tau.proc <-  1/(sig.proc*sig.proc)
  sig.proc ~ dunif(0,10) #dnorm(0,1)T(0,)
  tau.obs <-  1/(sig.obs*sig.obs)
  sig.obs ~ dunif(0,10) #dnorm(0,1)T(0,)
  
  # Nested random effects site/cliff
  for(i in 1:ncliffs) { eps.cliff[i] ~ dnorm( eps.site[site[i]],tau.cliff) } 
  for(j in 1:nsites) { eps.site[j]~ dnorm(0, tau.site)}
  tau.cliff  <- 1/(sig.cliff*sig.cliff) 
  sig.cliff ~ dunif(0,10) #dnorm(0,1)T(0,)
  tau.site  <- 1/(sig.site*sig.site) 
  sig.site ~ dunif(0,10) #dnorm(0,1)T(0,) 
  
  #### Data model ###### 
  for(k in 1 : ncounts) { 	
    lN.est[k] <- lmu.N[yr[k],state[k]] + eps.cliff[cliff[k]] 
    logtheta[k] ~ dnorm (lN.est[k], tau.obs)
    log(theta[k]) <- logtheta[k]
    count[k] ~ dpois(theta[k]) 
  }  # k
  for(s in 1 : nstates) { 
    # Baseline year and strata expectations ######
    lmn1st[s] <- log(mn1st[s]) # add a small constant to prevent log(zero) which will result in failure to run
    lmu.N[1,s] ~ dnorm(lmn1st[s], 0.1) # priors for 1st yr abundance are observed values
    # Dynamics
    for( y in startyear[s] : (endyear[s]-1)) { 
      lmu.N[y+1,s] <- lmu.N[y,s] + r[y,s] 
      r[y,s] ~ dnorm(mean.r, tau.proc)
      lambda[y,s] <- exp(r[y,s])
    } # y
  } # s
  
  # derived parameters
  rrate=mean(r[,2])
  mprate=mean(r[,1])
  ratediff=rrate-mprate
  r1=mean(r[1:3,2])
  r2=mean(r[6:9,2])
  rdiff=r2-r1
  # 
  mp1=mean(r[1:3,1])
  mp2=mean(r[6:9,1])
  mpdiff=mp2-mp1
  sdiff1=r1-mp1
  sdiff2=r2-mp2
  
  # Population sizes on real scale
  for (y in startyear[1] : endyear[1]) {
    MP.est[y] <- exp(lmu.N[y, 1])
    }
  for (y in startyear[2] : endyear[2]) {
    R.est[y] <- exp(lmu.N[y, 2])
  } 
} 
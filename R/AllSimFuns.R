Simpara <- function(nsites, 
                    p.prop, 
                    nreps,
                    Beta.min,
                    Beta.max,
                    addCovar = FALSE,
                    t.phi_mean, 
                    t.phi_sd,
                    sim.alpha = "uniform",
                    sdltheta = 0.5,
                    Ntheta = 1,
                    seed = 123){
  ### simulate parameters (mu, phi and theta) in NB based on design
  ### flag
  set.seed(seed = seed)
  flag = rep(FALSE, nsites)
  idx.TP = sample(1:nsites, size = round(nsites*p.prop))
  flag[idx.TP] = TRUE
  
  #### phi 
  set.seed(seed = seed)
  tmp = exp(rnorm(3*nsites, mean = t.phi_mean*(-5), sd = t.phi_sd*(1.5)))
  Phi0 = sample(tmp[tmp < 1], nsites, replace = FALSE)
  Phi0[Phi0 > 1] = 0.9999###?1
  
  # #### theta
  ################ modified on Mar 16, 2021
  ### 1. One theta: Same theta for IP, Input cross all conditions
  if(Ntheta == 1){
    set.seed(seed = seed)
    Theta0.uni = exp(rnorm(nsites,
                           #mean = 6.32 + 0.96*log(Phi0),
                           mean = 4.414 + 0.865*log(Phi0),
                           sd = sdltheta))
    Theta0.CtrlInput = Theta0.TrtInput = 
      Theta0.CtrlIP = Theta0.TrtIP = Theta0.uni
  }else if(Ntheta == 2){
    # ### 2. Two theta: Same theta_1 for Input, same theta_2 for IP cross all conditions
    set.seed(seed = seed)
    Theta0.Input = exp(rnorm(nsites,
                             mean = (4.414 + 0.865*log(Phi0))/2,
                             sd = sdltheta))
    set.seed(seed = seed)
    Theta0.IP = exp(rnorm(nsites,
                          mean = 4.414 + 0.865*log(Phi0),
                          sd = sdltheta))
    
    Theta0.CtrlInput = Theta0.TrtInput = Theta0.Input
    Theta0.CtrlIP = Theta0.TrtIP = Theta0.IP
  }else if(Ntheta == 3){
    # ### 3. Three theta: Same theta for IP and input in Ctrl, but 
    # ###               but different theta for IP and input in Trt
    set.seed(seed = seed)
    Theta0.Ctrl = exp(rnorm(nsites,
                            mean = 4.414 + 0.865*log(Phi0),
                            sd = sdltheta))
    set.seed(seed = seed)
    Theta0.TrtInput = exp(rnorm(nsites,
                                mean = (4.414 + 0.865*log(Phi0))*0.6,
                                sd = sdltheta))
    set.seed(seed = seed)
    Theta0.TrtIP = exp(rnorm(nsites,
                             mean = (4.414 + 0.865*log(Phi0))*0.8,
                             sd = sdltheta))
    
    Theta0.CtrlInput = Theta0.CtrlIP = Theta0.Ctrl
    Theta0.TrtInput = Theta0.TrtInput
    Theta0.TrtIP = Theta0.TrtIP
    
  }else if(Ntheta == 4){
    #### 4. Fourth theta: All phi are different
    set.seed(seed = seed)
    Theta0.CtrlInput = exp(rnorm(nsites,
                                 mean = (4.414 + 0.865*log(Phi0))*1.5,
                                 sd = sdltheta))
    set.seed(seed = seed)
    Theta0.CtrlIP = exp(rnorm(nsites,
                              mean = (4.414 + 0.865*log(Phi0))*2.25,
                              sd = sdltheta))
    set.seed(seed = seed)
    Theta0.TrtInput = exp(rnorm(nsites,
                                mean = (4.414 + 0.865*log(Phi0)),
                                sd = sdltheta))
    set.seed(seed = seed)
    Theta0.TrtIP = exp(rnorm(nsites,
                             mean = (4.414 + 0.865*log(Phi0))*1.5,
                             sd = sdltheta))
  }
  
  #################### combine all theta
  Theta0 = list(Input = cbind( matrix(rep(Theta0.CtrlInput, nreps[1]), 
                                      ncol = nreps[1], byrow = FALSE), 
                               matrix(rep(Theta0.TrtInput, nreps[2]), 
                                      ncol = nreps[2], byrow = FALSE)),
                IP = cbind( matrix(rep(Theta0.CtrlIP, nreps[1]), 
                                   ncol = nreps[1], byrow = FALSE), 
                            matrix(rep(Theta0.TrtIP, nreps[2]), 
                                   ncol = nreps[2], byrow = FALSE)))
  colnames(Theta0$Input) = c(paste0("Ctrl_input_rep", 1:nreps[1]),
                             paste0("Trt_input_rep", 1:nreps[2]))
  colnames(Theta0$IP) = c(paste0("Ctrl_IP_rep", 1:nreps[1]),
                          paste0("Trt_IP_rep", 1:nreps[2]))
  # #################
  
  
  
  #### mu
  if(sim.alpha == "uniform"){
    set.seed(seed = seed)
    alpha0 = runif(nsites, 0, 1) #runif(nsites, 0, 1) #rnorm(nsites, 1, 0.5) #rnorm(nsites, 0.5, 0.05) #runif(nsites, 0, 1)
  } else if(sim.alpha == "norm"){
    set.seed(seed = seed)
    #alpha0 = rnorm(nsites, 1.5, 0.4) 
    alpha0 = rnorm(nsites, 1.42, 0.5) 
    
  }else if(sim.alpha == "real"){
    #load("./data/Baseline_logitMu.rda")
    #alpha0 = sample(logit.Mu_pooled[, 1], size = nsites, replace = FALSE)
    #alpha0 = sample(logit.Mu_unpooled[[1]], size = nsites, replace = FALSE)
    load("./data/cbctx2wk_AlphaBeta.rda")
    set.seed(seed = seed)
    alpha0 = sample(alpha, size = nsites, replace = FALSE)
  }
  
  
  set.seed(seed = seed)
  beta0 = rep(0, nsites)
  beta0[idx.TP[1:(nsites*p.prop/2)]] = runif(nsites*p.prop/2, min = Beta.min, 
                                             max = Beta.max) # 1
  beta0[idx.TP[(nsites*p.prop/2 + 1):(nsites*p.prop)]] = runif(nsites*p.prop/2,
                                                               min = -Beta.max, max = -Beta.min) # -1
  
  if(addCovar){
    set.seed(seed)
    #gam0 = runif(nsites, 0.1, 0.5) 
    #### modified on April 2, 2021
   # gam0 = runif(nsites, 0.5, 1)  # old
    gam0 = runif(nsites, 1, 1.5)  # new
    ####
    
    set.seed(seed)
    age = runif(sum(nreps), 0, 1)
    design = data.frame(predictor = rep(c("Ctrl", "Trt"), nreps), 
                        Age = age)
    model = ~1 + predictor + Age
    model.matrix(model, design)
    
    eta0 = t(model.matrix(model, design)%*% t(cbind(alpha0, beta0, gam0))) ### methylation level in both groups
    mu0 = exp(eta0)/(1 + exp(eta0))
  }else{
    design = data.frame(predictor = rep(c("Ctrl", "Trt"), nreps))
    model = ~1 + predictor
    model.matrix(model, design)
    eta0 = t(model.matrix(model, design)%*% t(cbind(alpha0, beta0))) ### methylation level in both groups
    mu0 = exp(eta0)/(1 + exp(eta0))
  }
  
  colnames(mu0) = c(paste0("Ctrl_", 1:nreps[1]), paste0("Trt_", 1:nreps[2]))
  
  print(mu0[1777, ])
  
  if(addCovar){
    return(list(mu = mu0, phi = Phi0, theta = Theta0, 
                model = model, design = design, flag = flag, 
                alpha = alpha0, beta = beta0, gam = gam0))
  }else{
    return(list(mu = mu0, phi = Phi0, theta = Theta0, 
                model = model, design = design, flag = flag, 
                alpha = alpha0, beta = beta0))
  }
  
}

SimDat <- function(nreps, 
                   nsites, 
                   mu, 
                   phi,
                   theta,
                   logN_lambda = FALSE,
                   seed = 123){
  ### simulate data based on mu, phi and theta
  #### nreps: a vector of length G, containing the number of replicate in all G groups 
  #### mu: nsite*sum(nreps) matrix, denoting the methylation level in all samples, 
  ####     the order of samples are group 1, group2, group3,...
  #### phi: a vector of length nsites, denoting the dispersion of methylation for all sites
  #### theta: a list containing scale parameter of gamma for all sites, cross all conditions
  
  
  set.seed(seed = seed)
  ##### use mu to generate lambda
  call.rgamma <- function(x, n){
    #rgamma(n, shape = x[1:n], scale = x[length(x)])
    rgamma(n, shape = x[1:n], scale = x[(n+1):(2*n)])
  }
  
  plx = cbind((1- mu)*(1/phi - 1), theta$Input)
  lambda_x = t(apply(plx, 1, call.rgamma, n = sum(nreps)))
  
  ply = cbind(mu*(1/phi - 1), theta$IP)
  lambda_y = t(apply(ply, 1, call.rgamma, n = sum(nreps)))
  
  print((lambda_y/(lambda_x + lambda_y))[1777, ])

 
  ########### used for sensitivity analysis
  if(logN_lambda){
    call.norm <- function(x,n){
      exp(rnorm(n, mean = x[1:n], sd = x[length(x)]))
    }
     #par(mfrow = c(2,2))
     #hist(log(lambda_x), 100, xlim = c(0,8), main = "lambda.x ~ Gamma_x")
     #hist(log(lambda_y), 30, xlim = c(0,8), main = "lambda.y ~ Gamma_y")
     
    set.seed(seed = seed)
    ### randomly sample lambda from log-normal distribution
    # mu.x = log(lambda_x); sig.x = sqrt(rowVars(lambda_x))
    # tmp.x = cbind(mu.x, sig.x)
    ml.x = matrix(NA, nrow = nrow(lambda_x), ncol = sum(nreps))
    sdl.x = rep(NA, nrow(lambda_x))
    for (i in 1:nrow(lambda_x)) {
      ml.x[i,] = log(lambda_x[i,] + 0.00001) - (1/2)*log(var(lambda_x[i, ])/(mean(lambda_x[i,]))^2 + 1)
      sdl.x[i] = sqrt(log(var(lambda_x[i, ])/(mean(lambda_x[i,]))^2 + 1))
    }
    tmp.x = cbind(ml.x, sdl.x)
    lambda_x = t(apply(tmp.x, 1, call.norm, n = sum(nreps)))
    
    ###
    # mu.y = log(lambda_y); sig.y = sqrt(rowVars(lambda_y))
    # tmp.y = cbind(mu.y, sig.y)
    set.seed(seed = seed)
    ml.y = matrix(NA, nrow = nrow(lambda_y), ncol = sum(nreps))
    sdl.y = rep(NA, nrow(lambda_y))
    for (i in 1:nrow(lambda_y)) {
      ml.y[i,] = log(lambda_y[i,] + 0.00001) - (1/2)*log(var(lambda_y[i, ])/(mean(lambda_y[i,]))^2 + 1)
      sdl.y[i] = sqrt(log(var(lambda_y[i, ])/(mean(lambda_y[i,]))^2 + 1))
    }
    tmp.y = cbind(ml.y, sdl.y)
    lambda_y = t(apply(tmp.y, 1, call.norm, n = sum(nreps)))
    
    print((lambda_y/(lambda_x + lambda_y))[1777, ])
    
   # hist(log(lambda_x), 100, xlim = c(0,8), main = "lambda.x ~ logN_x")
   # hist(log(lambda_y), 30, xlim = c(0,8), main = "lambda.y ~ logN_y")
  }
  ###########
  
  ### size factor
  sf.x = runif(sum(nreps), 0.5, 1)
  sf.y = runif(sum(nreps), 0.5, 1)
  ### sample IP and control based on lambda_x, and lambda_y
  para.x = sweep(lambda_x, 2, sf.x, FUN = "*")
  para.y = sweep(lambda_y, 2, sf.y, FUN = "*")
  ### input ~ poison(para.x), 
  ### IP ~ poison(para.y)
  count.Input = matrix(rpois(n = nsites*sum(nreps), as.vector(para.x)), 
                       nrow = nsites, byrow = FALSE)
  count.IP = matrix(rpois(n = nsites*sum(nreps), as.vector(para.y)),
                    nrow = nsites, byrow = FALSE)
  
  tmp = sweep(count.IP, 2, sf.y, FUN = "/")/(sweep(count.IP, 2, sf.y, FUN = "/")+
                                               sweep(count.Input, 2, sf.x, FUN = "/"))
  #print(c(mean(tmp[1777, 1:nreps[1]]), mean(tmp[1777, (nreps[1] + 1):sum(nreps)])))
  print(tmp[1777, ])
  counts = matrix(0, nrow = nsites, ncol = 2*sum(nreps))
  ix = seq(1, ncol(counts), by = 2)
  counts[, ix] = count.Input
  ix = seq(2, ncol(counts), by = 2)
  counts[, ix] = count.IP
  num.na = sum(is.na(counts))
  if(length(num.na) > 0)
    counts[is.na(counts)] = rpois(num.na, 1e+06)
  
  colnames(counts) = 1:ncol(counts)
  colnames(counts)[seq(1, ncol(counts), 2)] = c(paste0("Ctrl_Input_Rep_",  1:nreps[1]),
                                                paste0("Trt_Input_Rep_",  1:nreps[2]))
  colnames(counts)[seq(2, ncol(counts), 2)] = c(paste0("Ctrl_IP_Rep_",  1:nreps[1]),
                                                paste0("Trt_IP_Rep_",  1:nreps[2]))
  sf = rep(0, 2*sum(nreps))
  sf[seq(1, ncol(counts), 2)] = sf.x
  sf[seq(2, ncol(counts), 2)] = sf.y
  return(list(counts = counts, sf = sf, 
              lambda_y = lambda_y, lambda_x = lambda_x,
              para.x = para.x, para.y = para.y))
}



######################################
Simpara.common <- function(nsites, 
                           p.prop, 
                           nreps,
                           Beta.min,
                           Beta.max,
                           addCovar = FALSE,
                           t.phi_mean, 
                           t.phi_sd,
                           sim.alpha = "uniform",
                           seed = 123){

  ### flag
  set.seed(seed = seed)
  flag = rep(FALSE, nsites)
  idx.TP = sample(1:nsites, size = round(nsites*p.prop))
  flag[idx.TP] = TRUE
  
  #### phi 
  set.seed(seed = seed)
  tmp = exp(rnorm(3*nsites, mean = t.phi_mean*(-5), sd = t.phi_sd*(1.5)))
  Phi0 = sample(tmp[tmp < 1], nsites, replace = FALSE)
  Phi0[Phi0 > 1] = 0.9999
  
  
  ### mu0
  if(sim.alpha == "uniform"){
    set.seed(seed = seed)
    alpha0 = runif(nsites, 0, 1) 
  } else if(sim.alpha == "norm"){
    set.seed(seed = seed)
    alpha0 = rnorm(nsites, 1.5, 0.4) 

  }else if(sim.alpha == "real"){
    # load("./data/cbctx2wk_AlphaBeta.rda")
    # set.seed(seed = seed)
    # alpha0 = sample(alpha, size = nsites, replace = FALSE)
    load("./data/simPara_GSE114150.rda")
    alpha0 = alpha
  }
  
  
  set.seed(seed = seed)
  beta0 = rep(0, nsites)
  beta0[idx.TP[1:(nsites*p.prop/2)]] = runif(nsites*p.prop/2, min = Beta.min, 
                                             max = Beta.max) # 1
  beta0[idx.TP[(nsites*p.prop/2 + 1):(nsites*p.prop)]] = runif(nsites*p.prop/2,
                                                               min = -Beta.max, max = -Beta.min) # -1
  
  if(addCovar){
    set.seed(seed)
    gam0 = runif(nsites, 1, 1.5)  

    
    set.seed(seed)
    age = runif(sum(nreps), 0, 1)
    design = data.frame(predictor = rep(c("Ctrl", "Trt"), nreps), 
                        Age = age)
    model = ~1 + predictor + Age
    model.matrix(model, design)
    
    eta0 = t(model.matrix(model, design)%*% t(cbind(alpha0, beta0, gam0))) ### methylation level in both groups
    mu0 = exp(eta0)/(1 + exp(eta0))
  }else{
    design = data.frame(predictor = rep(c("Ctrl", "Trt"), nreps))
    model = ~1 + predictor
    model.matrix(model, design)
    eta0 = t(model.matrix(model, design)%*% t(cbind(alpha0, beta0))) ### methylation level in both groups
    mu0 = exp(eta0)/(1 + exp(eta0))
  }
  
  colnames(mu0) = c(paste0("Ctrl_", 1:nreps[1]), paste0("Trt_", 1:nreps[2]))
  

  if(addCovar){
    return(list(mu = mu0, phi = Phi0,
                model = model, design = design, flag = flag, 
                alpha = alpha0, beta = beta0, gam = gam0))
  }else{
    return(list(mu = mu0, phi = Phi0, 
                model = model, design = design, flag = flag, 
                alpha = alpha0, beta = beta0))
  }
  
  
}

simDat.common <- function(nreps, 
                          nsites, 
                          mu, 
                          phi,
                          seed = 123){
  ### simulate data based on parameters
  set.seed(seed)
  load("./data/cb6wkvs2wk_Candidates.rda")
  X.real = Candidates$Counts[, seq(1, ncol(Candidates$Counts), 2)]
  idx = sample(1:nrow(X.real), nsites)
  lambda_x = sweep(X.real[idx,], 2, Candidates$sf[seq(1, length(Candidates$sf), 2)], FUN = "/")
  
  ### re-sample X
  set.seed(seed)
  sf.x = runif(sum(nreps), 0.5, 1)
  para.x = sweep(lambda_x, 2, sf.x, FUN = "*")
  count.Input = matrix(rpois(n = nsites*sum(nreps), as.vector(para.x)),
                       nrow = nsites, byrow = FALSE)
  
  
  # count.Input = X.real[idx, ]; sf.x = Candidates$sf[seq(1, length(Candidates$sf), 2)]
  # lambda_x = sweep(count.Input, 2, sf.x, FUN = "/")
  # 

  
  ### simulate p
  set.seed(seed = seed)
  call.rbeta <- function(x, n){
    rbeta(n, shape1 = x[1:n],  shape2 = x[(n+1):(2*n)])
    }
  
  beta.para = cbind(mu*(1/phi - 1), (1- mu)*(1/phi - 1))
  P = t(apply(beta.para, 1, call.rbeta, n = sum(nreps)))
  P[P > 0.99] = 0.99
  lambda_y = lambda_x * P/(1-P)
  
  
  
  ### size factor from real data
  set.seed(seed)
  sf.y = runif(sum(nreps), 0.5, 1) #Candidates$sf[seq(2,length(Candidates$sf), 2)] #1?why no real data here
  para.y = sweep(lambda_y, 2, sf.y, FUN = "*")
  count.IP = matrix(rpois(n = nsites*sum(nreps), as.vector(para.y)),
                    nrow = nsites, byrow = FALSE)
  
  # IP.norm = sweep(count.IP, 2, sf.y, FUN = "/")
  # Input.norm = sweep(count.Input, 2,sf.x, FUN = "/")

  
  ### write X and Y into one count matrix
  counts = matrix(0, nrow = nsites, ncol = 2*sum(nreps))
  ix = seq(1, ncol(counts), by = 2)
  counts[, ix] = count.Input
  ix = seq(2, ncol(counts), by = 2)
  counts[, ix] = count.IP
  colnames(counts) = 1:ncol(counts)
  colnames(counts)[seq(1, ncol(counts), 2)] = c(paste0("Ctrl_Input_Rep_",  1:nreps[1]),
                                                paste0("Trt_Input_Rep_",  1:nreps[2]))
  colnames(counts)[seq(2, ncol(counts), 2)] = c(paste0("Ctrl_IP_Rep_",  1:nreps[1]),
                                                paste0("Trt_IP_Rep_",  1:nreps[2]))
  
  
  sf = rep(0, 2*sum(nreps))
  sf[seq(1, ncol(counts), 2)] = sf.x
  sf[seq(2, ncol(counts), 2)] = sf.y
  
  return( list(counts = counts, sf = sf, P = P,
               lambda_y = lambda_y, lambda_x = lambda_x))
}





simDat.common2 <- function(nreps, 
                          nsites, 
                          mu, 
                          phi,
                          seed = 123){
  ### simulate data based on parameters
  load("./data/simPara_GSE114150.rda")

  set.seed(seed)
  disper_x = exp(rnorm(nsites, mean = -2.3, sd = 0.8))
  
  set.seed(seed)
  sf.x = runif(sum(nreps), 0.5, 1)
  sf.y = runif(sum(nreps), 0.5, 1)
  para.x = sweep(mu_x, 2, sf.x, FUN = "*")#2? mu_x is lambda?
  count.Input = matrix(rnbinom(n = nsites*sum(nreps), mu = as.vector(para.x), 
                               size = rep(1/disper_x, sum(nreps))), 
                       nrow = nsites, byrow = FALSE)
  
  
  
  ### simulate p
  call.rbeta <- function(x, n){
    rbeta(n, shape1 = x[1:n],  shape2 = x[(n+1):(2*n)])
  }
  set.seed(seed = seed)
  beta.para = cbind(mu*(1/phi - 1), (1- mu)*(1/phi - 1))
  P = t(apply(beta.para, 1, call.rbeta, n = sum(nreps)))
  P[P > 0.99] = 0.99
  mu_y = mu_x * P/(1-P)
  
  
  set.seed(seed = seed)
  para.y = sweep(mu_y, 2, sf.y, FUN = "*")
 # disper_y = exp(-0.36 + 0.325*log(disper_x))
  count.IP = matrix(rnbinom(n = nsites*sum(nreps), mu = as.vector(para.y), 
                     size = rep(1/disper_x, sum(nreps))), 
             nrow = nsites, byrow = FALSE)
  
  # 
  #  IP.norm = sweep(count.IP, 2, sf.y, FUN = "/")
  #  Input.norm = sweep(count.Input, 2,sf.x, FUN = "/")
  # Ratio =  (IP.norm/(IP.norm+Input.norm))
  #  
  
  ### write X and Y into one count matrix
  counts = matrix(0, nrow = nsites, ncol = 2*sum(nreps))
  ix = seq(1, ncol(counts), by = 2)
  counts[, ix] = count.Input
  ix = seq(2, ncol(counts), by = 2)
  counts[, ix] = count.IP
  colnames(counts) = 1:ncol(counts)
  colnames(counts)[seq(1, ncol(counts), 2)] = c(paste0("Ctrl_Input_Rep_",  1:nreps[1]),
                                                paste0("Trt_Input_Rep_",  1:nreps[2]))
  colnames(counts)[seq(2, ncol(counts), 2)] = c(paste0("Ctrl_IP_Rep_",  1:nreps[1]),
                                                paste0("Trt_IP_Rep_",  1:nreps[2]))
  
  
  sf = rep(0, 2*sum(nreps))
  sf[seq(1, ncol(counts), 2)] = sf.x
  sf[seq(2, ncol(counts), 2)] = sf.y
  
  return( list(counts = counts, sf = sf, P = P,
               mu_y = mu_y, mu_x = mu_x, 
               disper_x = disper_x, disper_y = disper_y))
}




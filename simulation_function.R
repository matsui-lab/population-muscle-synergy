

simEMG_envelope_multi  <- function(ref,nsim=100,
                                   n_jump = 2,
                                   amp_rate = 1,
                                   noise_rate = 1,
                                   plot.it=FALSE){
  
  basis <- splines::bs(x = seq_along(ref),knots = seq(1,length(ref),by = 30))
  bslm <- lm(ref ~ 0 + basis)
  cf <- bslm$coefficients
  #recov <- basis %*% cf
  prob <- rbeta(length(cf),1,2)
  jump_idx <- sample(length(cf),n_jump,prob = prob)
  
  sim_coef <- cf
  mu <- rep(0,length(sim_coef))
  if(n_jump > 0){
    jump_mu <- sample(c(1,-1),n_jump,replace = TRUE)*rnorm(n_jump,amp_rate*max(abs(ref)),1)
    mu[jump_idx] <- jump_mu
    #sim_coef[jump_idx] <- sim_coef[jump_idx] + rnorm(length(jump_idx),jump_mu,sd(ref)*noise_rate)
    sim_coef <- sim_coef + rnorm(length(sim_coef),mu,sd(ref)*noise_rate)
  }
  
  
  sim_coef2 <- matrix(rep(sim_coef,nsim),nrow=length(sim_coef),ncol=nsim)
  jump_mu <- sample(c(1,-1),n_jump,replace = TRUE)*rnorm(n_jump,amp_rate*max(abs(ref)),1)
  
  for(i in 1:ncol(sim_coef2)){
    #sim_coef[jump_idx,i] <- sim_coef[jump_idx,i] + rnorm(length(jump_idx),jump_mu,sd(ref)*noise_rate)
    sim_coef2[,i] <- sim_coef2[,i] + rnorm(nrow(sim_coef2),0,sd(ref)*noise_rate)
  }
  
  simdat <- basis %*% sim_coef2
  
  if(plot.it){
    matplot(simdat,type="l")
    matplot(rowMeans(simdat),type="l",add=TRUE,lwd=3,col=1)
    matplot(ref,type="l",add=TRUE,lwd=3)
  }
  simdat
}

simEMG_envelope_multi_simple  <- function(x,nsim=100,noise_rate=1,plot.it=FALSE){
  
  basis <- splines::bs(x = seq_along(x),knots = seq(1,length(x),by = 30))
  bslm <- lm(x ~ 0 + basis)
  cf <- bslm$coefficients
  
  mu <- rep(0,length(cf))
  sim_coef <- matrix(0,nrow=length(cf),ncol=nsim)
  for(i in 1:ncol(sim_coef)){
    sim_coef[,i] <- cf + rnorm(nrow(sim_coef),mu,sd(x)*noise_rate)
  }
  simdat <- basis %*% sim_coef
  simdat
}



simEMG_envelope  <- function(ref,
                             n_jump = 2,
                             amp_rate = 1,
                             noise_rate = 1,
                             plot.it=FALSE){
  
  basis <- splines::bs(x = seq_along(ref),knots = seq(1,length(ref),by = 30))
  bslm <- lm(ref ~ 0 + basis)
  cf <- bslm$coefficients
  #recov <- basis %*% cf
  prob <- rbeta(length(cf),1,2)
  jump_idx <- sample(length(cf),n_jump,prob = prob)
  
  sim_coef <- cf
  mu <- rep(0,length(sim_coef))
  if(n_jump > 0){
    jump_mu <- sample(c(1,-1),n_jump,replace = TRUE)*rnorm(n_jump,amp_rate*max(abs(ref)),1)
    mu[jump_idx] <- jump_mu
    #sim_coef[jump_idx] <- sim_coef[jump_idx] + rnorm(length(jump_idx),jump_mu,sd(ref)*noise_rate)
    sim_coef <- sim_coef + rnorm(length(sim_coef),mu,sd(ref)*noise_rate)
  }
  
  
  # sim_coef <- matrix(rep(cf,nsim),nrow=length(cf),ncol=nsim)
  # jump_mu <- sample(c(1,-1),n_jump,replace = TRUE)*rnorm(n_jump,amp_rate*max(abs(ref)),1)
  # 
  # for(i in 1:ncol(sim_coef)){
  #   sim_coef[jump_idx,i] <- sim_coef[jump_idx,i] + rnorm(length(jump_idx),jump_mu,sd(ref)*noise_rate)
  # }
  
  simdat <- basis %*% sim_coef
  
  if(plot.it){
    matplot(simdat,type="l")
    matplot(rowMeans(simdat),type="l",add=TRUE,lwd=3,col=1)
    matplot(ref,type="l",add=TRUE,lwd=3)
  }
  simdat
}

sim_amp <- function(ref,
                    n_jump = 4,
                    n_burst = 100,
                    noise_rate = 2,
                    plot.it=FALSE){
  
  len <- length(ref)
  signal0 <- signal <- ref * sign(rnorm(len,0,1)) * rgamma(len,10,10)
  if(n_jump == 0){
    return(signal0)
  }
  
  prob <- rbeta(len,1,2)
  jump_idx <- sample(len,n_jump,prob = prob)
  jump_idx_st <- jump_idx - ceiling(n_burst/2)
  jump_idx_ed <- jump_idx + ceiling(n_burst/2)
  jump_idx_st[jump_idx_st <= 0] <- 1
  jump_idx_ed[jump_idx_ed >= len] <- len
  #  print(jump_idx_st)
  #  print(jump_idx_ed)
  
  for(ii in seq_along(jump_idx_st)){
    st <- jump_idx_st[ii]
    ed <- jump_idx_ed[ii]
    Sd <- noise_rate*sd(signal0)
    eps <- rnorm(length(st:ed),0,Sd)
    signal[st:ed] <- signal[st:ed] + eps
    #    print(signal[st:ed])
    #    print(paste("is.na.noise",any(is.na(eps))))
  }
  #  print(length(signal))
  
  if(plot.it){
    matplot(signal,type="l",col=2)
    matplot(signal0,type="l",add=TRUE)
  }
  signal
}



sim_amp_multi <- function(ref,nsim=100,
                          n_jump = 4,
                          n_burst = 100,
                          noise_rate = 2,
                          plot.it=FALSE){
  
  len <- length(ref)
  signal0 <- signal <- ref * sign(rnorm(len,0,1)) * rgamma(len,10,10)
  signal0 <- signal <- sapply(1:nsim,function(ii)ref * sign(rnorm(len,0,1)) * rgamma(len,10,10))
  if(n_jump == 0){
    return(signal0)
  }
  
  prob <- rbeta(len,1,2)
  jump_idx <- sample(len,n_jump,prob = prob)
  jump_idx_st <- jump_idx - ceiling(n_burst/2)
  jump_idx_ed <- jump_idx + ceiling(n_burst/2)
  jump_idx_st[jump_idx_st <= 0] <- 1
  jump_idx_ed[jump_idx_ed >= len] <- len
  #  print(jump_idx_st)
  #  print(jump_idx_ed)
  
  for(ii in seq_along(jump_idx_st)){
    for(jj in 1:nsim){
      st <- jump_idx_st[ii]
      ed <- jump_idx_ed[ii]
      Sd <- noise_rate*sd(signal0)
      eps <- rnorm(length(st:ed),0,Sd)
      signal[st:ed,jj] <- signal[st:ed,jj] + eps
      #    print(signal[st:ed])
      #    print(paste("is.na.noise",any(is.na(eps))))
    }
  }
  #  print(length(signal))
  
  if(plot.it){
    matplot(signal,type="l",col=2)
    matplot(signal0,type="l",add=TRUE)
  }
  signal
}



sim_sEMG <- function(refmat,
                     n_jump_envelope = 2,
                     amp_rate = 1,
                     noise_rate_envelope = 0.5,
                     noise_rate_amp = 1.5,
                     n_jump_amp = 4,
                     n_burst = 1,
                     artifact_rate = 0.3,
                     outlier_rate = 0.3,
                     n_trial = 100){
  
  # Subject patterns
  
  n_artifact <- ceiling(n_trial*artifact_rate)
  n_normal <- n_trial - n_artifact
  
  # Trial patterns
  
  n_normal_outlier <- ceiling(outlier_rate*n_normal)
  n_normal <- n_normal - n_normal_outlier
  
  n_artifact_outlier <- ceiling(outlier_rate*n_artifact)
  n_artifact <- n_artifact - n_artifact_outlier
  
  sim_signal <- array(NA,dim = c(nT,ncol(refmat),n_trial))
  
  sim_envelope1 <- apply(refmat,2,
                         function(xx)simEMG_envelope(xx,
                                                     n_jump = 0,
                                                     amp_rate = amp_rate,
                                                     noise_rate = noise_rate_envelope))
  sim_envelope2 <- apply(refmat,2,
                         function(xx)simEMG_envelope(xx,
                                                     n_jump = n_jump_envelope,
                                                     amp_rate = amp_rate,
                                                     noise_rate = noise_rate_envelope))
  
  
  for(i in 1:n_normal){
    sim_signal[,,i] <- apply(sim_envelope1,2,
                             function(xx)sim_amp(xx,
                                                 n_jump = 0,
                                                 noise_rate = noise_rate_amp))
  }
  
  st <- n_normal+1
  ed <- n_normal + n_normal_outlier
  
  for(i in st:ed){
    sim_signal[,,i] <- apply(sim_envelope1,2,
                             function(xx)sim_amp(xx,
                                                 n_jump = n_jump_amp,
                                                 n_burst = n_burst,
                                                 noise_rate = noise_rate_amp))
  }
  
  st <- n_normal + n_normal_outlier + 1
  ed <- n_artifact + n_normal + n_normal_outlier
  
  for(i in st:ed){
    
    sim_signal[,,i] <- apply(sim_envelope2,2,
                             function(xx)sim_amp(xx,
                                                 n_jump = 0,
                                                 noise_rate = noise_rate_amp))
  }
  
  st <- n_artifact + n_normal + n_normal_outlier + 1
  ed <- n_artifact_outlier + n_artifact + n_normal + n_normal_outlier
  
  for(i in st:ed){
    sim_signal[,,i] <- apply(sim_envelope2,2,
                             function(xx)sim_amp(xx,
                                                 n_jump = n_jump_amp,
                                                 n_burst = n_burst,
                                                 noise_rate = noise_rate_amp))
  }
  
  
  sim_signal
}


t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               maxColorValue = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}




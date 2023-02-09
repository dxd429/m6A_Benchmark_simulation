Plot_ROC <- function(flag, Pvals, models,ltypes, plabel, cols = NULL, leg, axis_size = 2.5, label_size = 3.5,
                     mgp = c(5,2,0), lwd = 6, cex = 4, tcl = -0.8, lwd_x = 4, lwd_y = 4){
  library(ROCR)

  ROC.DE <- function(DE.gs, pval) {
    pred = prediction(1-pval, DE.gs)
    perf = performance(pred,"tpr","fpr")
    perf
  }
  nmodel = ncol(Pvals)
  if(is.null(cols)){
    cols = 1:nmodel
  }
  for(i in 1:nmodel){
    idx = is.na( Pvals[, i])
    roc = ROC.DE(flag[!idx], Pvals[!idx, i])
    symb<- rep(NA, length(roc@x.values[[1]]))
    symb[seq(1, length(roc@x.values[[1]]), by = round(length(roc@x.values[[1]])/10))]<- plabel[i]
    symb[1]<- NA
    if (i == 1){
      xlim=c(0,1)
      par(mgp=mgp)
      
      plot(roc, axes = FALSE, xlab="1-Specificity", ylab="Sensitivity", type = "o", xlim = xlim, lty= ltypes[i], col = alpha(cols[i],1), lwd = lwd, pch = symb,yaxis.at=c(0,0.2,0.4,0.6,0.8,1),
           xaxis.cex.axis=axis_size,yaxis.cex.axis=axis_size, tcl = tcl, cex.lab = label_size, cex = cex,yaxis.lwd=lwd_y, xaxis.lwd = lwd_x)
      

    }else if(i > 1){
      par(mgp=mgp)

      plot(roc, axes = FALSE, xlab="1-Specificity", ylab="Sensitivity",xaxt = "n", type = "o", add = TRUE, col = alpha(cols[i],1), lty = ltypes[i], lwd = lwd, pch = symb,yaxis.at=c(0,0.2,0.4,0.6,0.8,1),
           xaxis.cex.axis=axis_size,yaxis.cex.axis=axis_size, tcl = tcl, cex.lab = label_size, cex = cex,yaxis.lwd=lwd_y, xaxis.lwd = lwd_x)
     
    }
  }
  if(leg)
    legend("bottomright", legend= models, col = alpha(cols,1), lty=ltypes, pch = plabel, cex = 1.8)
  axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
       labels = paste0(c('', '', '', '', '', "")),
       cex.axis = axis_size, lwd = lwd_x, tck = tcl) 
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
       labels = paste0(c('', '', '', '', '', "")),
       cex.axis = axis_size,lwd = lwd_y, tck = tcl) 
}



Plot_PreRecall <- function(flag, Pvals, models, cols = NULL, 
                           ltypes = 1:6, leg){
  library(ROCR)
  PrecRec.DE <- function(DE.gs, pval) {
    pred = prediction(1-pval, DE.gs)
   # perf = performance(pred,"tpr","fpr")
    perf = performance(pred,"prec","rec")
    perf
  }
  nmodel = ncol(Pvals)
  if(is.null(cols)){
    cols = 1:nmodel
  }
  for(i in 1:nmodel){
    idx = is.na( Pvals[, i])
    PrecRec = PrecRec.DE(flag[!idx], Pvals[!idx, i])
    
    if (i == 1){
      xlim=c(0.05,1)
      plot(PrecRec, xlim = xlim, ylim=c(0,1), lty = ltypes[i], 
           col = cols[i])
    }else if(i > 1){
      plot(PrecRec, add = TRUE, col = cols[i],
           lty = ltypes[i], cex.axis = 1.1)
    }
  }
  if(leg)
    legend("bottomleft", legend= models, col = cols, lty=ltypes, cex = 1)
}


MLE.NB.phi <- function(phi, y){
  #### find the MLE of NB's dispersion 
  s = 0
  for (j in 1:length(y)) {
    tmp = y[j]/phi - (1/phi + y[j])*(mu[j] /(1+mu[j]*phi)) +
      (1/phi^2)*log(1+mu[j]*phi) - sum(1/(phi + seq(0,(y[j]-1))*phi^2))
    s = s+tmp
  }
  
  return(s)
}


# log.NB <- function(t, y){
#   log.lik = sum(dnbinom(y, size = exp(-t) + 1, mu = y, log = TRUE))
#   log.lik
# }



Plot_PrecTop <- function(flag, Pvals, Ntop = 1000,
                         models, cols = NULL, ltypes = 1:6,
                         yylim,
                         leg){
  ### calculate the proportion of true DM sites among the top ranked sites
  Ranks = seq(100, Ntop, 100)
  cum.prop = matrix(NA, nrow = length(Ranks), ncol = ncol(Pvals))
  for (i in 1:ncol(Pvals)) {
    thispval = Pvals[, i]
    tmp.flag = flag[order(thispval)]
    ### claculate cumulative proportion 
    cum.prop[, i] = cumsum(tmp.flag)[Ranks]/ Ranks
  }
  
  matplot(cum.prop, type = "l", xaxt = "n", 
          xlab = "Top ranked regions", ylab = "True DM m6A regions",
          col = cols, lty = ltypes, lwd = 2.2,
          cex.axis = 1.1,
          ylim = yylim)
  if(leg){
    legend("bottomleft", legend = models, cex = 0.6, col= cols,  lty = ltypes)
  }
  
  
  axis(1, at = seq(length(Ranks)/5, length(Ranks), by = length(Ranks)/5),
       labels = paste0(seq(length(Ranks)/5, length(Ranks), by = length(Ranks)/5), "00"),
       cex.axis = 1.1) 
 
}



Plot_PrecTop1 <- function(flag, Pvals, Ntop = 1000,
                          models, cols = NULL, ltypes,
                          yylim,
                          leg){
  ### calculate the proportion of true DM sites among the top ranked sites
  Ranks = seq(100, Ntop, 100)
  cum.prop = matrix(NA, nrow = length(Ranks), ncol = length(Pvals))
  for (m in 1:length(Pvals)) {
    thisPval = as.matrix(Pvals[[m]])
    this.cum.prop = matrix(NA, nrow = length(Ranks), ncol = ncol(thisPval))
    for (isim in 1:ncol(thisPval)) {
      pval.sim = thisPval[, isim]
      tmp.flag = flag[order(pval.sim)]
      ### claculate cumulative proportion 
      this.cum.prop[, isim] = cumsum(tmp.flag)[Ranks]/ Ranks
    }
    if(ncol(this.cum.prop) > 1){
      cum.prop[, m] = rowMeans(this.cum.prop, na.rm = TRUE)
    }else{
      cum.prop[, m] = this.cum.prop
    }
  }
  
  
  matplot(cum.prop, type = "l", xaxt = "n", 
          xlab = "Top ranked regions", ylab = "True DM m6A regions",
          col = cols, lty = ltypes, lwd = 2.2,
          cex.axis = 1.1,
          ylim = yylim)
  if(leg){
    legend("bottomleft", legend = models, cex = 0.9, col= cols,  lty = ltypes,
           bty = "n")
  }
  
  
  axis(1, at = seq(length(Ranks)/5, length(Ranks), by = length(Ranks)/5),
       labels = paste0(seq(length(Ranks)/5, length(Ranks), by = length(Ranks)/5), "00"),
       cex.axis = 1.1) 
  
}



#FDR_t(t0 = threshold0, pval = PVAL.QNB[input_idx[[j]][[it]]], flag = flag[input_idx[[j]][[it]]])

### FDR and type I error under threshold t0
FDR_t <- function(t0, pval, fdr = NULL, flag){
  if(length(fdr) == 0)
    fdr = p.adjust(pval, method = "fdr")
  R = sum(fdr < t0, na.rm = TRUE)
  V = sum(fdr < t0 & (!flag), na.rm = TRUE)
  return(V/R)
}

TypeIEr_t <- function(p0, pval, flag){
  V = sum(pval < p0 & (!flag), na.rm = TRUE)
  return(V/sum(!flag))
}




TransPval <- function(stats){
  ##### 1. mixture normal
  res.EM = mix.2norm.onlysd(Y = stats[!is.na(stats)], pi = 0.1)
  sd0.mix = res.EM$sd1
  pval.mix = 2*(1- pnorm(abs(stats), mean = 0, sd = sd0.mix) )
  fdr.mix = p.adjust(pval.mix, method = "fdr")
  #sum(pval.mix < 1e-3 , na.rm = TRUE);sum(fdr.mix<0.05, na.rm = TRUE)
  
  #### 2. truncated normal
  sd0.TrunN = Uniroot.truncNsd(Y = stats, a = -2, b = 2)
  pval.TrunN = 2*(1- pnorm(abs(stats), mean = 0, sd = sd0.TrunN) )
  fdr.TrunN = p.adjust(pval.TrunN, method = "fdr")
  #sum(pval.TrunN < 1e-3, na.rm = TRUE);sum(fdr.TrunN<0.05, na.rm = TRUE)
  
  
  #### 3. local fdr 
  # library(locfdr)
  # idx = (!is.na(stats))
  # res.locfdr = locfdr(zz = stats[idx], nulltype = 1,plot = 0)
  # Locfdr[idx] = res.locfdr$fdr
  # emp.locfdr = sum(Locfdr < 0.2 &!flag, na.rm = TRUE)/sum(Locfdr<0.2, na.rm = TRUE)
  library(locfdr)
  Global.fdr = rep(NA, length(stats))
  # idx = (!is.na(stats))
  # normstat = stats[idx]
  # #normstat = normstat[normstat > -60]
  # normstat[normstat < -20] = -20
  # normstat[normstat > 20] = 20
  # fdrres = locfdr(normstat, plot = 0)
  # fdr.global = numeric(length(normstat))
  # xx = fdrres$mat[, "x"]
  # leftbreaks = xx - (xx[2] - xx[1])/2
  # leftbreaks[1] = min(normstat)
  # rightbreaks = xx + (xx[2] - xx[1])/2
  # rightbreaks[length(rightbreaks)] = max(normstat)
  # for (ii in 1:length(normstat)) {
  #   ind = ((leftbreaks <= (-1) * abs(normstat[ii])) | 
  #            (rightbreaks >= abs(normstat[ii])))
  #   F1l = sum(fdrres$mat[ind, "p1f1"])
  #   Fl = sum(fdrres$mat[ind, "f"])
  #   fdr.global[ii] = 1 - F1l/Fl
  # }
  # Global.fdr[idx] = fdr.global
  # sum(Global.fdr < 0.05, na.rm = TRUE)
  
  Infe = list(pval.mix = pval.mix, fdr.mix = fdr.mix, sd0.mix = sd0.mix,
                    pval.TrunN = pval.TrunN, fdr.TrunN = fdr.TrunN, sd0.TrunN = sd0.TrunN,
                    Locfdr = Global.fdr)
  
}




#### permFDR
permFDR <- function(realStat,  permStat){
  ### calculate FDR for each site based on permutation
  real.ord = order(abs(realStat), decreasing = TRUE)
  n = length(realStat)
  perm.FDR = perm.fdr = rep(NA, n)
  for (ir in 1:length(real.ord)) {
    if(!is.na(realStat[real.ord[ir]])){
      thisstat = abs(realStat[real.ord[ir]])
      n.large = sum(abs(permStat) >= thisstat, na.rm = TRUE)
      perm.fdr[ir] = min(n.large/ir, 1)
    }
  }
  
  perm.FDR[real.ord] = perm.fdr
  return(perm.FDR)
}

logit <- function(x){
  x[x==0] = 0.0001
  x[x==1] = 0.999
  log(x/(1-x))
}


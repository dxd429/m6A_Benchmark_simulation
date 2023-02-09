##################################################  fun for DESeq2
DiffPeak.DESeq2 <- function(counts, nreps, 
                            model,
                            design,
                            sf = NULL) {
  library(DESeq2)
  # design = data.frame(Reps = c(rep(paste0("Rep", 1:nreps[1]),each = 2),
  #                              rep(paste0("Rep", 1:nreps[2]),each = 2)),
  #                     IP = rep(c("Input", "IP"), sum(nreps)), 
  #                     Trt = rep(c("Ctrl", "Trt"), 2*nreps ))
  # model = ~ Reps + IP + Trt + IP*Trt
  # model.matrix(model, design)
  # 
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = design,
                                design = model)
  
  if(length(sf) ==  0){
    ## compute size factor from all counts
    sf = colSums(counts)
    sf = sf/median(sf)
  }
  counts.norm = sweep(counts, 2, sf, FUN = "/")
  sizeFactors(dds) = sf
  dds <- DESeq(dds, test = "Wald")
  
  ###
  #dds <- nbinomWaldTest(dds)
                           
  ###
  
  resultsNames(dds) # name of each variable in the design matrix used for DESeq
  res = DESeq2::results(dds, name = "IPIP.TrtTrt") 
  dat = data.frame(pval = res$pvalue, stat = res$stat)
  return(dat)
}





#########################################################################################  funs for FET-HMM

# high spatial resolution differential analysis of RNA methylome with 
# fisher exact test and Hidden Markov Model
rhtestHMM <- function(untreated_ip, untreated_input, treated_ip, treated_input, 
                      untreated_ip_total, untreated_input_total, treated_ip_total, treated_input_total,
                      Pi=matrix(c(0.9,0.1,0.1,0.9),byrow=TRUE, nrow=2),
                      delta=c(0.5,0.5),
                      pm=list(prob=c(0.1, 0.9)),
                      threshold_fdr=0.05,
                      gene_label=rep(1,length(untreated_ip)),
                      minimal_count_fdr =10,mode="DIRECT"){
  # parameters
  # untreated_ip_total: an integer, sequencing depth of untreated ip sample   
  # untreated_input_total: ... untreated input sample                         
  # treated_ip_total: ... treated ip sample                                   
  # treated_input_total: ... treated input sample                             
  
  # untreated_ip:  a vector of integers, vector length equal to the number of binding sites, each element contains the number of reads within a binding site for the untreated ip sample
  # untreated_input:  ... of the untreated input sample
  # treated_ip:  ... of the treated ip sample
  # treated_input:  ... of the treated input sample
  
  # Pi:initial state transition probability matrix of HMM. default: matrix(c(0.9,0.1,0.1,0.9),byrow=TRUE, nrow=2)
  # delta:initial state probability of HMM. default: c(0.5,0.5)
  # pm:initial emission probability of HMM. default: c(0.1,0.9), i.e., the probability of observing a differential methylation state on an undifferential site is 0.1; and the probability for observing a differential methylation state on an differential methylation site is 0.9.
  # threshold_fdr:the threshold used to convert bltest states to observations of HMM model. default: 0.05, i.e., if a site is differential methylated with significance level 0.05, than it is considered a differential methylation loci in HMM model.
  # gene_label:if the counts are from multiple genes, then HMM is built on each gene. default: rep(1,length(untreated_ip)), i.e., the counts are all from the same gene.
  
  # minimal_count_fdr:an integer threshold, only the loci with reads more than this number are subjected for fdr calculation. default: 10
  # mode:a string,which specifies the mode of using HMM default:default:"DIRECT",which means to apply Estep of EM algorithm after estimating transition matrix and initial probability for all genes.
  # Alternative:"EM", which means to smooth the significance level(p-value) directly; "BEL",which means to use the Bernoulli process as the observation process of HMM.
  # required library
  require(exomePeak)
  require(HiddenMarkov)
  # rhtest analysis
  result<-rhtest(untreated_ip, untreated_input, treated_ip, treated_input, untreated_ip_total,
                 untreated_input_total, treated_ip_total, treated_input_total,
                 minimal_count_fdr = minimal_count_fdr)
  labelrhtest<- (result$log.fdr<log(threshold_fdr))
  if(mode=="DIRECT"){
    deltaestimate<-c(1-mean(labelrhtest),mean(labelrhtest))
    length<-length(treated_ip)
    state0<-c()
    p01<-c()
    p00<-c()
    state1<-c()
    p11<-c()
    p10<-c()
    # Pi
    state0 <- (labelrhtest[-length]==FALSE)
    p01<-mean(labelrhtest[-1][state0]) # 0 -> 1
    p00<-1- mean(labelrhtest[-1][state0]) # 0 -> 0
    state1 <- (labelrhtest[-length]==TRUE)
    p11<-mean(labelrhtest[-1][state1]) # 1 -> 1
    p10<-1 - mean(labelrhtest[-1][state1]) # 1 -> 0
    piestimate<-matrix(c(p00,p01,p10,p11), nrow = 2, ncol = 2, byrow = 2,dimnames = NULL)
    # HMM section
    pp<-exp(result$log.p)
    prob<-cbind(pp,1-pp)
    pen <- colSums(prob^2)
    prob[,1] <- prob[,1]/pen[1]
    prob[,2] <- prob[,2]/pen[2]
    unique_gl = unique(gene_label)
    pos <- cbind(untreated_ip*0 + 1,untreated_ip*0)
    #Estepresult<-list()
    for (i in 1:length(unique_gl)) {
      id = which(gene_label == unique_gl[i])
      m<-length(deltaestimate)
      len<-length(id)
      y<-.Estep_prob(piestimate, deltaestimate, m,len,prob[id,],fortran = TRUE, fwd.only = FALSE)
      pos[id,]<-y$u
    }
    # save HMM result
    m = untreated_ip + untreated_input + treated_ip + treated_input
    IDreads= which(m >minimal_count_fdr)
    #IDnot=which(m<= minimal_count_fdr)
    log.fdr<-c(rep(0,length(treated_ip)))
    log.p<-c(rep(0,length(treated_ip)))
    log.fdr[IDreads]=log(p.adjust(pos[IDreads,1], method = "fdr"))
    # remove infinity
    log.fdr=pmax(log.fdr,-1000)
    log.p=pmax(log(pos[,1]),-1000)
    #log.p[IDnot]=0
    # find peaks
    tmp <- ctest(IP=untreated_ip+treated_ip,INPUT=untreated_input+treated_input,
                 TOTAL_IP=untreated_ip_total+treated_ip_total,
                 TOTAL_INPUT=untreated_input_total+treated_input_total)
    tmp2 <- (tmp$log.fdr > log(0.05))
    # modify result
    log.fdr[tmp2] <- 0
    log.p[tmp2] <- 0
    result$log.fdr[tmp2]<-0
    result$log.p[tmp2]<-0
    fdr_fisher_sort <-sort(result$log.fdr)
    p_fisher_sort <-sort(result$log.p)
    log.fdr<- fdr_fisher_sort[rank(log.fdr)]
    log.p<- p_fisher_sort[rank(log.p)]
  }
  if(mode=="BEL"){
    # HMM section
    unique_gl = unique(gene_label)
    # initialization the posterior probability
    pos <- untreated_ip*0 + 1;
    for (i in 1:length(unique_gl)){
      id = which(gene_label == unique_gl[i])
      labelrhtest_i <- labelrhtest[id]
      pos[id] <- .HMM_single(labelrhtest_i,Pi,delta,pm)
    }
    # save HMM result
    log.fdr=log(p.adjust(pos, method = "fdr"))
    # adjust only sig reads count testing result
    m = untreated_ip + untreated_input + treated_ip + treated_input
    ID= which(m > minimal_count_fdr)
    log.fdr_sig=log(p.adjust(pos[ID], method = "fdr"))
    log.fdr[ID] = log.fdr_sig
    # remove infinity
    log.fdr=pmax(log.fdr,-1000)
    log.p=pmax(log(pos),-1000)
    
  }else if(mode=="EM"){
    # HMM section
    pp<-exp(result$log.p)
    prob<-cbind(pp,1-pp)
    # add the penalty to the p-values
    pen <- colSums(prob^2)
    prob[,1] <- prob[,1]/pen[1]
    prob[,2] <- prob[,2]/pen[2]
    unique_gl = unique(gene_label)
    # using the EM method(Baum-Welch Algorithm)
    rhtestHMM<-.Mstep_prob(Pi,delta,m,len, prob)
    # save HMM result
    log.fdr=log(p.adjust(rhtestHMM$gama, method = "fdr"))
    # adjust only sig reads count testing result
    m = untreated_ip + untreated_input + treated_ip + treated_input
    ID= which(m > minimal_count_fdr)
    log.fdr_sig=log(p.adjust(rhtestHMM$gama[ID], method = "fdr"))
    log.fdr[ID] = log.fdr_sig  
    # remove infinity
    log.fdr=pmax(log.fdr,-1000)
    log.p=pmax(log(rhtestHMM$gama),-1000)
  }
  # save result
  DIFF=list(log.fdr=log.fdr,log.p=log.p,log.fc=result$log.fc, rhtest_result = result)
  return(DIFF)
}
##########################################################################
# subfunction1
.HMM_single <- function(labelrhtest,Pi,delta,pm) {
  aa<-length(labelrhtest)
  pn <- list(size=rep(1,aa))
  x <- dthmm(labelrhtest, Pi, delta, "binom", pm, pn,discrete=TRUE)
  log <- capture.output({
    y <- BaumWelch(x);
  })
  pos <- y$u[,1]
  return(pos)
}
###############################################################################
#subfuction2
#define the Estep.prob function
.Estep_prob<-function(Pi, delta,m,len, pr, fortran = TRUE, fwd.only = FALSE){
  y <- forwardback.dthmm(Pi, delta, pr, fortran = TRUE, fwd.only = FALSE)
  logbeta <- y$logbeta
  logalpha <- y$logalpha
  LL <- y$LL
  u <- exp(logalpha + logbeta - LL)
  v <- array(NA, dim = c(len - 1, m, m))
  for (k in 1:m) {
    logPi <- matrix(log(Pi[, k]), byrow = TRUE, nrow = len -1, ncol = m)
    logPbeta <- matrix(log(pr[-1,k]) + logbeta[-1, k], byrow = FALSE, 
                       nrow = len - 1, ncol = m)
    v[, , k] <- logPi + logalpha[-len, ] + logPbeta - LL
  }
  v <- exp(v)
  return(list(u = u, v = v, LL = LL))
}
###############################################################################
#subfuction3(define the Mstep function)
.Mstep_prob<-function(pistart, deltastart,m,len, prob){
  pinew <- pistart
  deltanew <- deltastart
  iter<-1
  LL<-c()
  LL[1]=0
  # perform the E-step and M-step of Baum-Welch Algorithm
  repeat 
  { 
    m <- nrow(pinew)
    len<-length(untreated_ip)
    tsik<-c()
    tsij<-matrix(data = NA, nrow = 2, ncol = 2)
    y<-.Estep_prob(pinew, deltanew, m,len,prob, fortran = TRUE, fwd.only = FALSE)
    gama<- y$u
    ksi<- y$v
    LL[iter+1]<-y$LL
    # HMM M-step ,tisk & tsij are the expectation of times of transition.
    for (i in 1:2){
      tsik[i]=sum(gama[1:len-1,i])
      for (j in 1:2){
        tsij[i,j]=sum(ksi[1:len-1,i,j])}
    }
    # update delta and pi
    for (i in 1:2){
      deltanew[i]=gama[1,i]
      for (j in 1:2){
        pinew[i,j]=tsij[i,j]/tsik[i]}
    }
    print(iter)
    print(LL[iter])
    if(abs(LL[iter+1]-LL[iter])<(0.1^4))
    {
      break;
    }
    iter=iter+1;
  }
  #|| iter>100
  #the HMM converges
  iter<-iter
  gama<- y$u
  ksi<- y$v
  LLfinal<-y$LL
  return(list(gama = gama[,1], iter= iter, LLfinal = LLfinal))
}






########################################################################### simple linear model
DiffPeak.lm <- function(Ratio, design){
  ### Ratio: a matirx containing ratios between normalized IP and input in all groups
  ### design: design matrix
  Ratio[Ratio == 0] = 0.01
  Ratio[Ratio == 1] = 0.99
  Z = as.matrix(log(Ratio/(1-Ratio)))
  
  R = matrix(NA, nrow = nrow(Z), ncol = ncol(design))
  Cov.R = vector("list", length = nrow(Z))
  for (i in 1:nrow(Ratio)) {
    R[i, ] = (t(solve(t(design)%*%design)%*%t(design)%*%Z[i,])) ### coefficients, one row for one site
    tmp_i = sum((Z[i, ]%*%(diag(1, nrow(design)) - design%*%solve(t(design)%*%design)%*%t(design) ))*t(Z[i, ]))
    Cov.R[[i]] = tmp_i*solve(t(design)%*%design)
    colnames( Cov.R[[i]]) = 
      rownames( Cov.R[[i]]) = colnames(design)
  }
  
  R = as.data.frame(R)
  colnames(R) = colnames(design)
  
  #### calculate statistics
  se.trt = rep(NA, nrow(Ratio))
  for (isite in 1:nrow(Ratio)) {
    tmp = Cov.R[[isite]]
    se.trt[isite] = sqrt(tmp["TrtTrt", "TrtTrt"])
  }
  stat = R$TrtTrt/se.trt
  
  return(list(Coef = R, Cov = Cov.R, 
              stat = stat))
}


##################################################  fun for exomePeak

Bltest <- function(control_ip, treated_ip, control_input, treated_input){
  library(exomePeak)
  control_ip_total <- sum(colSums(control_ip))
  control_input_total <- sum(colSums(control_input))
  treated_ip_total <- sum(colSums(treated_ip))
  treated_input_total <- sum(colSums(treated_input))
  print("start test...")
  tmpResult <- do.call(cbind.data.frame, exomePeak::bltest(rowSums(control_ip), 
                                                           rowSums(control_input),
                                                           rowSums(treated_ip),
                                                           rowSums(treated_input),
                                                           control_ip_total,
                                                           control_input_total,
                                                           treated_ip_total,  
                                                           treated_input_total) ) 
  return( data.frame(logFC = tmpResult[,"log.fc"], pvalue = exp(tmpResult[,"log.p"]), 
                     fdr = exp(tmpResult$log.fdr) ) )
}



##################################################  fun for MeTDiff
MeTDiffTest <-  function( control_ip, treated_ip, control_input, treated_input ){
  testResult <- foreach(i = 1:nrow(control_ip) , .combine = rbind )%do%{
    x <- t( as.matrix(control_ip[i,]) )
    y <- t( as.matrix(control_input[i,]) )
    xx <- t( as.matrix(treated_ip[i,]) )
    yy <- t( as.matrix(treated_input[i,]) )
    xxx = cbind(x,xx)
    yyy = cbind(y,yy)
    
    logl1 <- MeTDiff:::.betabinomial.lh(x,y+1)
    logl2 <- MeTDiff:::.betabinomial.lh(xx,yy+1)
    logl3 <- MeTDiff:::.betabinomial.lh(xxx,yyy+1)
    tst <- (logl1$logl+logl2$logl-logl3$logl)*2
    pvalues <- 1 - pchisq(tst,2)
    log.fc <- log( (sum(xx)+1)/(1+sum(yy)) * (1+sum(y))/(1+sum(x)) ) 
    
    data.frame( logFC = log.fc, pvalue =  pvalues )
  }
  testResult$fdr <- p.adjust(testResult$pvalue, method = "fdr")
  return(testResult)
}


################################################## all funs from RADAR

DiffPeak.RADAR <- function(counts, variable){
  
  ### Counts: order of columns are: 
  #### Ctrl.input_1, ..., Ctrl.input_N.reps[1], Trt.input_1, ..., Trt.input_N.reps[2],
  #### Ctrl.ip_1, ..., Ctrl.ip_N.reps[1], Trt.ip_1, ..., Trt.ip_N.reps[2]
  
  ### variable: contain variables in the design, where the first column corresponds to the factor of interest, or the testing factor. 
  #### The rest columns are covariates
  library(RADAR)
  
  RADAR.dat <- MeRIP()
  
  ### added on Feb 2, 2021 in order to avoid error when transforming MeRIP class to MeRIP.RADAR
  ### this txdb will never be used in simulation
  TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene -> txdb
  RADAR.dat@geneModel = exonsBy(txdb,by="gene")
  ###
  
  
  RADAR.dat@reads = as.matrix(counts) 
  rownames(RADAR.dat@reads) = as.character(paste0("Peak_", 1:nrow(counts) ))
  RADAR.dat@samplenames = as.character(paste0("s", 1:nrow(variable)) )
  
  RADAR.dat <- normalizeLibrary( RADAR.dat, boxPlot = FALSE)
  if(length(levels(variable$predictor)) == 0){
    variable$predictor = as.factor(variable$predictor)
  }
  RADAR.dat <- adjustExprLevel( RADAR.dat )
  variable(RADAR.dat) <- variable
  RADAR.dat <- filterBins(RADAR.dat, minCountsCutOff = 15)
  RADAR.dat <- diffIP(RADAR.dat)
  library(readr)
  
  
  res = matrix(NA, nrow = nrow(counts), ncol = ncol(RADAR.dat@test.est))
  rm.idx = parse_number(setdiff(rownames(RADAR.dat@ip_adjExpr), 
                                rownames(RADAR.dat@ip_adjExpr_filtered)))
  if(length(rm.idx) > 0){
    res[-rm.idx, ] = RADAR.dat@test.est
  }else{
    res = RADAR.dat@test.est
  }
  colnames(res) = colnames(RADAR.dat@test.est) 
  rownames(res) = rownames(RADAR.dat@ip_adjExpr)
  res = as.data.frame(res)
  return(res)
  
  # pval.RADAR = rep(NA, nrow(counts))
  # if(length(rm.idx) > 0){
  #   pval.RADAR[-rm.idx] = RADAR.dat@test.est[, "p_value"]
  # }else{
  #   pval.RADAR = RADAR.dat@test.est[, "p_value"]
  # }
  # return(pval.RADAR)
}





normalizeLibrary.RADAR <- function(object, boxPlot = TRUE){
  ### this functin is directly copied from RADAR on github
  ## load data from input
  m6A <- object@reads[,(1+length(object@samplenames)):(2*length(object@samplenames))]
  input <- object@reads[,1:length(object@samplenames)]
  colnames(input) <- colnames(m6A) <-  object@samplenames
  object@geneBins <- geneBins<- geneBins(object)
  
  ## Get input geneSum (gene level quantification)
  geneSum <- NULL
  for(i in 1:ncol(input) ){
    y <- input[,i]
    gene.sum <- tapply(y,geneBins$gene,sum)
    geneSum <- cbind(geneSum,gene.sum)
  }
  colnames(geneSum) <- object@samplenames
  
  size.input <- DESeq2::estimateSizeFactorsForMatrix(geneSum)
  
  ## compute normalized input data
  norm.input <-t( t(input) / size.input )
  geneSum.norm <- t ( t(geneSum)/size.input)
  
  
  ## estimate enrichment using top IP count bins
  ave.ip <- rowMeans(m6A)
  ave.top <- order(ave.ip,decreasing = T)[1:round(0.01*length(ave.ip)[1])]
  
  ## Get the gene level input count for corresponding bins
  geneCounts.window <- geneSum.norm[geneBins[rownames(m6A),"gene"],]
  
  enrich <- as.data.frame(m6A[ave.top,]/geneCounts.window[ave.top,])
  enrich <- enrich[!apply(enrich,1, function(x){any(is.na(x)) | any(is.infinite(x))}),]
  
  size.enrich.deseq2 <- DESeq2::estimateSizeFactorsForMatrix(enrich[,1:length(object@samplenames)])
  
  ## calculate normzlied ip read count
  norm.ip <-t( t(m6A)/size.enrich.deseq2 )
  sizeFactor <- data.frame(input=size.input,ip=size.enrich.deseq2)
  
  object@geneSum <- geneSum.norm
  
  outObj <- as(object, "MeRIP.RADAR" )
  outObj@norm.ip <- norm.ip
  outObj@norm.input <- norm.input
  outObj@sizeFactor <- sizeFactor
  
  if(boxPlot){
    plot.new()
    par(mfrow=c(2,2))
    boxplot(log(geneSum[rowSums(geneSum)!=0,]+1),main = "INPUT")
    boxplot(log(geneSum.norm[rowSums(geneSum.norm)!=0,]+1),main = "Normalized INPUT")
    boxplot(log(enrich[rowSums(enrich)!=0,]+0.1), main = "IP (Estimated enrichment)")
    enrich.norm <- as.data.frame(norm.ip[ave.top,]/geneCounts.window[ave.top,])
    boxplot(log(enrich.norm[rowSums(enrich.norm)!=0,]+0.1), main = "Normalized IP (estimated enrichment)")
    par(mfrow=c(1,1))
  }
  
  return(outObj)
  
}





################################################## all funs from exomePeak2, essentially DESeq2
DiffPeak.exomePeak2 <- function(se,
                                txdb,
                                test_method,
                                p_cutoff,
                                exByGene,
                                bin_size,
                                alt_hypothesis,
                                lfc_threshold,
                                motif_based){
  #Set reference levels
  se$IP_input <- relevel(factor(se$IP_input),"input")
  se$Perturbation <- relevel(factor(se$Perturbation),"C")
  dds <- DESeqDataSet(se, ~ IP_input * Perturbation)
  
  normalizationFactors(dds) <- assays(se)[["sfm"]]
  
  #Fit differential models
  if(test_method == "DESeq2"){
    dds <- estimateDispersions(dds)
  }else{
    dispersions(dds) <- 0
  }
  
  dds <- nbinomWaldTest(dds)
  res <- results(dds,
                 altHypothesis = alt_hypothesis,
                 lfcThreshold = lfc_threshold)
  pvalue <- res$pvalue
  lfc <- res$log2FoldChange
  rm(dds)
  
  #Merge bins and annotate the resulting peaks
  pvalue[is.na(pvalue)] <- 1
  peak <- reducePeaks(rowRanges(se)[pvalue <= p_cutoff], exByGene)
  if(!motif_based) peak <- peak[sum(width(peak)) >= bin_size]
  exbg <- exonsBy(txdb, by = "gene")
  
  mcols(peak) <- makePeakAnnot(peak, se, res, exbg)
  return(peak)
}




################################################## all funs from DRME

DMEseq <- function(meth1,meth2,unmeth1,unmeth2) {
  # don't improve robustness, we don't add 1 to reads count
  s <- .sizeFactor2(cbind(meth1,meth2,unmeth1,unmeth2))
  s_t1 <- s[1:length(meth1[1,])]
  s_t2 <- s[(length(meth1[1,])+1):(length(cbind(meth1,meth2)[1,]))]
  s_c1 <- s[(length(cbind(meth1,meth2)[1,])+1):(length(cbind(meth1,meth2,unmeth1)[1,]))]
  s_c2 <- s[(length(cbind(meth1,meth2,unmeth1)[1,])+1):(length(cbind(meth1,meth2,unmeth1,unmeth2)[1,]))]
  
  # estimate probability of methylation under a condition
  p1 <- .estimateP(meth1, unmeth1, s_t1, s_c1)
  p2 <- .estimateP(meth2, unmeth2, s_t2, s_c2)
  p0 <- .estimateP(cbind(meth1,meth2), cbind(unmeth1,unmeth2), c(s_t1,s_t2), c(s_c1,s_c2))
  
  # estimate the abundance of feature
  q0 <- .estimateQ(cbind(meth1,meth2), cbind(unmeth1,unmeth2), c(s_t1,s_t2), c(s_c1,s_c2),p0)
  q1 <- .estimateQ(meth1,unmeth1,s_t1,s_c1,p0)
  q2 <- .estimateQ(meth2,unmeth2,s_t2,s_c2,p0)
  
  # estimate size e
  e1 <- q1/q0
  e2 <- q2/q0
  
  
  # calculate methylation reads count variance on common scale, condition 1
  w_t1 <-.calculateW(meth1,s_t1,e1)
  w_t2 <-.calculateW(meth2,s_t2,e2)
  w_c1 <-.calculateW(unmeth1,s_c1,e1)
  w_c2 <-.calculateW(unmeth2,s_c2,e2)
  
  # locfit 
  fit_t1 <- .locfitW(p1,q0,w_t1)
  fit_t2 <- .locfitW(p2,q0,w_t2)
  fit_c1 <- .locfitW(p1,q0,w_c1)
  fit_c2 <- .locfitW(p2,q0,w_c2)
  .plotDispersion(fit_t1,fit_t2)
  
  # calculate z
  z_t1 <- .calculateZ(q0,p1,s_t1,e1)
  z_t2 <- .calculateZ(q0,p2,s_t2,e2)
  z_c1 <- .calculateZ(q0,(1-p1),s_c1,e1)
  z_c2 <- .calculateZ(q0,(1-p2),s_c2,e2)
  
  # get estimate w
  w_fit_t1 <- .fittedW(p0,q0,fit_t1)
  w_fit_t2 <- .fittedW(p0,q0,fit_t2)
  w_fit_c1 <- .fittedW(p0,q0,fit_c1)
  w_fit_c2 <- .fittedW(p0,q0,fit_c1)
  
  # get estimate of upi
  ups_t1 <- pmax(w_fit_t1 - z_t1, 1e-8)
  ups_t2 <- pmax(w_fit_t2 - z_t2, 1e-8)
  ups_c1 <- pmax(w_fit_c1 - z_c1, 1e-8)
  ups_c2 <- pmax(w_fit_c2 - z_c2, 1e-8)
  
  # get all means
  mu_t1 <- (e1*q0*p0)%*%t(as.numeric(s_t1))
  mu_t2 <- (e2*q0*p0)%*%t(as.numeric(s_t2))
  mu_c1 <- (e1*q0*(1-p0))%*%t(as.numeric(s_c1))
  mu_c2 <- (e2*q0*(1-p0))%*%t(as.numeric(s_c2))
  
  # get all variance
  raw_t1 <- (e1%*%t(s_t1))^2*ups_t1
  raw_t2 <- (e2%*%t(s_t2))^2*ups_t2
  raw_c1 <- (e1%*%t(s_c1))^2*ups_c1
  raw_c2 <- (e2%*%t(s_c2))^2*ups_c2
  
  # put mu together
  mu1_t <- rowSums(mu_t1)
  mu2_t <- rowSums(mu_t2)
  mu1_c <- rowSums(mu_c1)
  mu2_c <- rowSums(mu_c2)
  
  # put size together
  size1_t <- (mu1_t^2)/rowSums(raw_t1)
  size2_t <- (mu2_t^2)/rowSums(raw_t2)
  size1_c <- (mu1_c^2)/rowSums(raw_c1)
  size2_c <- (mu2_c^2)/rowSums(raw_c2)
  
  # observation together
  t1 <- rowSums(meth1)
  t2 <- rowSums(meth2)
  c1 <- rowSums(unmeth1)
  c2 <- rowSums(unmeth2)
  t <- t1 + t2
  n1 <- t1 + c1
  n2 <- t2 + c2
  
  # go to test
  res <- .quadNBtest(t1,t,n1,n2,mu1_t,mu2_t,mu1_c,mu2_c,size1_t,size2_t,size1_c,size2_c)
  
  # add fc
  fc <- log(p1/p2)
  m1 <- rowSums(t(t(meth1)/s_t1))
  m2 <- rowSums(t(t(meth2)/s_t2))
  u1 <- rowSums(t(t(unmeth1)/s_c1))
  u2 <- rowSums(t(t(unmeth2)/s_c2))
  mfc <- log(m1)-log(m2)
  ufc <- log(u1)-log(u2)
  
  res2 <- data.frame(res,fc,mfc,ufc,q0)
  res <- res2[,c(1,3:7)]
  colnames(res) <- c("pvalue","mu","fc","meth.fc","unmeth.fc","q0")
  return(res)}

# Other Sub-functions
.sizeFactor <- function(n, useTotal=FALSE) {
  if (useTotal) {
    temp <- log(colSums(n))
    temp <- temp - mean(temp)
    s_size <- exp(temp)
  } else {
    n <- pmax(n,1e-5)
    log_n <- log(n)
    pseudo <- rowMeans(log_n)
    ratio <- log_n-pseudo
    s_size <- exp(apply(ratio,2,median)) }
  return(s_size)
}
.sizeFactor2 <- function(n) {
  temp <- log(colSums(n))
  temp <- temp - mean(temp)
  s_size <- exp(temp)
  return(s_size)
}
.estimateP <- function(meth, unmeth, size_t, size_c) {
  temp_t <- t(t(meth)/size_t)
  temp_c <- t(t(unmeth)/size_c)
  temp_n <- temp_t+temp_c
  p <- rowSums(temp_t)/rowSums(temp_n)
  p[is.na(p)] <- 0.5
  p[is.infinite(p)] <- 0.5
  return(p)
  # which(is.na(p))
  # which(is.infinite(p))  
}
.estimateQ <- function(meth,unmeth,size_t,size_c,p,useAll=FALSE){
  if (useAll) {
    temp_t <- t(t(meth)/size_t)
    temp_c <- t(t(unmeth)/size_c)
    temp_n <- temp_t+temp_c
    q <- rowSums(temp_n)/length(size_t)
  } else {
    temp_c <- t(t(unmeth)/size_c)
    qc <- rowMeans(temp_c)
    q <- qc/(1-p)
  }
  
  q[is.na(q)] <- 0
  return(q)
  # which(is.na(q))
  # which(is.infinite(q))  
}
.calculateW <- function(meth,size_t,e_t){
  temp_t <- t(t(meth)/size_t)
  q <- temp_t/e_t
  #  q[is.na(q)] <- 0
  resi <- q-rowMeans(q)
  w <- rowSums(resi^2)/(length(size_t)-1)
  w <- pmax(w,1e-8)
  return(w)
}
.locfitW <- function(p,q,w) {
  l <- log(q+1)
  data=data.frame(cbind(p,l,w))
  ID <- which(rowSums(is.na(data))>0)
  #  data <- data[-ID,]
  fit=locfit(w~lp(p,l),data=data,family="gamma")
  return(fit)
}
.calculateZ <- function(q,p,size,e){
  temp <- p*q/length(size)
  
  norow <- length(q)
  nocol <- length(size)
  temp2 <- matrix(1,nrow=norow, ncol=nocol )
  temp3 <- t(t(temp2)/size)/e
  
  z <- rowSums(temp3)*temp
  
  #  z[is.infinite(z)] <- 0
  #  z[is.na(z)] <- 0
  return(z)
}
.fittedW <- function(p,q,fit){ 
  library(locfit)
  l <- log(q+1)
  data=data.frame(cbind(p,l))
  w_fit <- predict(fit,data)
  return(w_fit)
}
.quadNBtest <- function(t1,t,n1,n2,mu1_t,mu2_t,mu1_c,mu2_c,size1_t,size2_t,size1_c,size2_c){
  nrows <- length(t)
  pval4 <- rep(1,nrows)
  pval2 <- rep(1,nrows)
  
  t2 <- t-t1
  for (irow in 1:nrows) {
    print(irow)
    trip <- t[irow]
    
    if (trip<1) {p <- NA} else {
      
      trip_t1 <- 0:t[irow]
      trip_t2 <- t[irow] - trip_t1
      trip_c1 <- n1[irow] - trip_t1
      trip_c2 <- n2[irow] - trip_t2
      
      
      p1 <- dnbinom(x=trip_t1, size=size1_t[irow], mu=mu1_t[irow], log = TRUE)
      p2 <- dnbinom(x=trip_t2, size=size2_t[irow], mu=mu2_t[irow], log = TRUE)
      p3 <- dnbinom(x=trip_c1, size=size1_c[irow], mu=mu1_c[irow], log = TRUE)
      p4 <- dnbinom(x=trip_c2, size=size2_c[irow], mu=mu2_c[irow], log = TRUE)
      
      #p4
      p <- p1+p2+p3+p4
      p <- p-max(p)
      p <- exp(p)/sum(exp(p))
      p <- min(1,2*min(sum(p[1:(t1[irow]+1)]),sum(p[(t1[irow]+1):(t[irow]+1)])))- p[(t1[irow]+1)]/2  
      pval4[irow] <- p
      
      #p2
      p <- p1+p2
      p <- p-max(p)
      p <- exp(p)/sum(exp(p))
      p <- min(1,2*min(sum(p[1:(t1[irow]+1)]),sum(p[(t1[irow]+1):(t[irow]+1)])))- p[(t1[irow]+1)]/2  
      pval2[irow] <- p
    }
    
  }
  
  mu <- (mu1_t+mu2_t+mu1_c+mu2_c)/4
  
  res <- data.frame(pval2,pval4,mu)
  return(res)
}
.plotDispersion <-function(fit_t1,fit_t2) {
  p <- rep(seq(from=0,to=1,by=0.01),10)
  q <- rep(seq(from=0,to=10,by=0.1),10)
  p <- matrix(p,nrow = 101, ncol = 101,byrow=TRUE)
  q <- matrix(q,nrow = 101, ncol = 101,byrow=FALSE)
  
  p <- p[1:10201]
  q <- q[1:10201]
  
  w1 <- .fittedW(p,q,fit_t1)
  w1 <- matrix(log(w1),nrow = 101, ncol = 101,byrow=FALSE)
  w2 <- .fittedW(p,q,fit_t2)
  w2 <- matrix(log(w2),nrow = 101, ncol = 101,byrow=FALSE)
  
  pdf("dispersion.pdf",height=4,width=7)
  par(mfrow=c(1,2))
  contour(seq(from=0,to=1,by=0.01),seq(from=0,to=10,by=0.1),w1)
  contour(seq(from=0,to=1,by=0.01),seq(from=0,to=10,by=0.1),w2)
  dev.off()
}





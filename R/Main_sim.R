##########################################################################
#### Main simulation workflow
##########################################################################
sim_ind<- 1 # Iteration Index
library(tictoc)
library(Matrix)
library(QNB)
library(exomePeak)
library(DESeq2)
library(RADAR)
source("simDat.real.R")
source("OtherDiffTools.R")
source("Utils_simulation.R")

library(miceadds)
source.all("./TRESS_funs/")

# N.sites = 10000 
# P.prop = 0.1 #c(0.05, 0.1, 0.2) #c(0.05, 0.1, 0.2)
# BETA = cbind(c(0.5, 1.5)) #cbind(c(0.3, 0.6), c(0.5, 1))
AddCov = FALSE
# T.phi_mean = 1#c(1, 0.8)
# T.phi_sd = 0.6 
N.reps.Ctrl = c(2,3,5,7,10)#c(2, 5) 
N.reps.Trt  = c(2,3,5,7,10)#c(2, 5) 
nsim = 1

PVAL_noCov = PVAL_Cov = vector("list",length = length(N.reps.Ctrl))
#paras_noCov = paras_Cov = vector("list",length = length(N.reps.Ctrl))
        for (i in 1:length(N.reps.Ctrl)) {
            N.reps = c(N.reps.Ctrl[i], N.reps.Trt[i])
            ###### start to simulate data
            res.sim = simDat.real(AddCov = AddCov,
                               nreps = N.reps,
                               seed = sim_ind)

            # res.para = Simpara.common(nsites = nsites,
            #                           nreps = N.reps,
            #                           p.prop = p.prop,
            #                           Beta.min = Beta.min,
            #                           Beta.max = Beta.max,
            #                           addCovar = AddCov,
            #                           t.phi_mean = t.phi_mean,
            #                           t.phi_sd = t.phi_sd,
            #                           sim.alpha = "real")
            flag = res.sim$flag
            model.matrix(res.sim$model, res.sim$design)
            stat.TRESS = pval.TRESS = pval.exome2 = pval.fethmm = pval.drme = 
              pval.qnb = pval.exome = pval.metdiff = pval.RADAR = matrix(NA, nrow = length(flag), ncol = nsim)
            Contrast = c(c(0, 1), rep(0, ncol(res.sim$design)-1))
 
           save(res.sim, file = paste0("/mnt/pan/SOM_PQHS_HXF155/daoyu/migrate/m6a/Evaluation/Evaluation/AltSim/res_sim_",N.reps.Ctrl[i], "_", N.reps.Trt[i], "_seed", sim_ind, ".rda"))
              # res.sim = simDat.common2(nreps = N.reps,
              #                          nsites = nsites,
              #                          mu = res.para$mu,
              #                          phi = res.para$phi,
              #                          seed = 1)
              
              #  Ratio = meRatio(res.sim$counts, res.sim$sf)
              #  Ratio[which(res.para$flag)[1:5],]
              #  res.sim$P[which(res.para$flag)[1:5], ]
              #  
              # phi.mom = rowVars(Ratio)/(rowMeans(Ratio)*(1-rowMeans(Ratio)))
              # # (res.sim$IP.norm/(res.sim$IP.norm+res.sim$Input.norm))[which(res.para$flag)[1:5],]
              # ##########################################################
              
              
              #### TRESS
              tic("TRESS")
              TRESS.DMR = CallDMRs.paramEsti(counts = res.sim$counts,
                                             sf = res.sim$sf,
                                             variable = res.sim$design,
                                             model = res.sim$model)
              TRESS.test = TRESS_DMRtest(DMR = TRESS.DMR, contrast = Contrast, nullModel = "standN")
             # stat.TRESS[, 1] = TRESS.test$stat
              pval.TRESS[, 1] = TRESS.test$pvalue
              bench_time<- toc()
              time_TRESS<- diff(as.numeric(bench_time[1:2]))
              
              # # ####
              # plot(res.sim$alpha, TRESS.DMR$Coef[,1], cex = 0.5, pch = 16)
              # abline(0,1,col = "red")
              # plot(res.sim$beta, TRESS.DMR$Coef[,2], cex = 0.5, pch = 16)
              # abline(0,1,col = "red")
              # plot(log(res.sim$phi), log(TRESS.DMR$shrkPhi),
              #      xlim = c(-8, 0), ylim= c(-8,0),cex = 0.3 )
              # abline(0,1,col = "red")
              # hist(log(res.sim$phi), xlim = c(-8, 0))
              # hist(log(TRESS.DMR$shrkPhi), xlim = c(-8, 0))
              # # ####
              
              
              
              #### exomePeak2, essentially desea2
              tic("exome2")
              design.exome2 = data.frame(Reps = c(rep(paste0("Rep", 1:N.reps[1]),each = 2),
                                                  rep(paste0("Rep", 1:N.reps[2]),each = 2)),
                                         IP = rep(c("Input", "IP"), sum(N.reps)),
                                         Trt = rep(c("Ctrl", "Trt"), 2*N.reps ))
              model.exome2 = ~ Reps + IP + Trt + IP*Trt
              model.matrix(model.exome2, design.exome2)
              res.exome2 = DiffPeak.DESeq2(counts = res.sim$counts, nreps = N.reps, 
                                           model = model.exome2, design = design.exome2) 
              pval.exome2[, 1] = res.exome2$pval
              bench_time<- toc()
              time_exome2<- diff(as.numeric(bench_time[1:2]))
              
              
              #### QNB
              tic("QNB")
              Ctrl_Input.id = grepl("Ctrl_Input", colnames(res.sim$counts))
              Ctrl_IP.id = grepl("Ctrl_IP", colnames(res.sim$counts))
              Trt_Input.id = grepl("Trt_Input", colnames(res.sim$counts))
              Trt_IP.id = grepl("Trt_IP", colnames(res.sim$counts))
              
              res.QNB = qnbtest(control_input = res.sim$counts[, Ctrl_Input.id],
                                control_ip = res.sim$counts[, Ctrl_IP.id],
                                treated_input = res.sim$counts[, Trt_Input.id],
                                treated_ip = res.sim$counts[, Trt_IP.id],
                                size.factor = list(control_input = res.sim$sf[Ctrl_Input.id ],
                                                   control_ip = res.sim$sf[Ctrl_IP.id],
                                                   treated_input = res.sim$sf[Trt_Input.id],
                                                   treated_ip = res.sim$sf[Trt_IP.id]
                                ),
                                mode = "per-condition",
                                output.dir = "./Results/QNB1")
              pval.qnb[, 1] = res.QNB$pvalue
              bench_time<- toc()
              time_qnb<- diff(as.numeric(bench_time[1:2]))
              
              
              #### exomePeak
              tic("exome")
              res.exomePeak = Bltest(control_ip = res.sim$counts[, Ctrl_IP.id] ,
                                     treated_ip = res.sim$counts[, Trt_IP.id],
                                     control_input = res.sim$counts[, Ctrl_Input.id],
                                     treated_input = res.sim$counts[, Trt_Input.id] )
              
              pval.exome[, 1] = res.exomePeak$pvalue
              bench_time<- toc()
              time_exome<- diff(as.numeric(bench_time[1:2]))
              
              #### MeTDiff
              tic("metdiff")
              res.Metdiff <- MeTDiffTest(control_ip = res.sim$counts[, Ctrl_IP.id] ,
                                         treated_ip = res.sim$counts[, Trt_IP.id],
                                         control_input = res.sim$counts[, Ctrl_Input.id],
                                         treated_input = res.sim$counts[, Trt_Input.id])
              pval.metdiff[, 1] = res.Metdiff$pvalue
              bench_time<- toc()
              time_metdiff<- diff(as.numeric(bench_time[1:2]))
              
              #### FET-HMM
              tic("fethmm")
              res.fethmm <- rhtestHMM(untreated_ip = rowSums(res.sim$counts[, Ctrl_IP.id]) ,
                                      untreated_input = rowSums(res.sim$counts[, Ctrl_Input.id]), 
                                      treated_ip = rowSums(res.sim$counts[, Trt_IP.id]), 
                                      treated_input = rowSums(res.sim$counts[, Trt_Input.id]), 
                                      untreated_ip_total = sum(rowSums(res.sim$counts[, Ctrl_IP.id])), 
                                      untreated_input_total = sum(rowSums(res.sim$counts[, Ctrl_Input.id])), 
                                      treated_ip_total = sum(rowSums(res.sim$counts[, Trt_IP.id])), 
                                      treated_input_total = sum(rowSums(res.sim$counts[, Trt_Input.id])),
                                      mode="DIRECT")
              pval.fethmm[, 1] = exp(res.fethmm$rhtest_result$log.p)
              bench_time<- toc()
              time_fethmm<- diff(as.numeric(bench_time[1:2]))
              
              
              #### DRME
              tic("drme")
              res.drme = DMEseq(unmeth1 = res.sim$counts[, Ctrl_Input.id],  
                                meth1 = res.sim$counts[, Ctrl_IP.id],
                                unmeth2 = res.sim$counts[, Trt_Input.id],
                                meth2  = res.sim$counts[, Trt_IP.id])
              pval.drme[, 1] = res.drme$pvalue
              bench_time<- toc()
              time_drme<- diff(as.numeric(bench_time[1:2]))
              
              #### RADAR
              tic("RADAR")
              radar.idx = c(seq(1, ncol(res.sim$counts), 2), seq(2, ncol(res.sim$counts), 2))
              radar.counts = res.sim$counts[, radar.idx]
              radar.var <- res.sim$design
              res.RADAR = DiffPeak.RADAR(counts = radar.counts, variable = radar.var)
              pval.RADAR[, 1] = res.RADAR$p_value
              bench_time<- toc()
              time_RADAR<- diff(as.numeric(bench_time[1:2]))
            
            
            #### save PVALS
            id = i
            if(AddCov){
              PVAL_Cov[[id]] = list(flag = flag,
                                    QNB = pval.qnb, 
                                    Time_QNB = time_qnb,
                                    exomePeak = pval.exome,
                                    Time_exomePeak = time_exome, 
                                    exomePeak2 = pval.exome2,
                                    Time_exomePeak2 = time_exome2,
                                    MeTDiff = pval.metdiff,
                                    Time_MeTDiff = time_metdiff,
                                    RADAR = pval.RADAR, 
                                    Time_RADAR = time_RADAR,
                                    FETHMM = pval.fethmm,
                                    Time_FETHMM = time_fethmm,
                                    DRME = pval.drme,
                                    Time_DRME = time_drme,
                                    TRESS = pval.TRESS,
                                    Time_TRESS = time_TRESS
                                    )
              
            }else{
              PVAL_noCov[[id]] = list(flag = flag,
                                    QNB = pval.qnb, 
                                    Time_QNB = time_qnb,
                                    exomePeak = pval.exome,
                                    Time_exomePeak = time_exome, 
                                    exomePeak2 = pval.exome2,
                                    Time_exomePeak2 = time_exome2,
                                    MeTDiff = pval.metdiff,
                                    Time_MeTDiff = time_metdiff,
                                    RADAR = pval.RADAR, 
                                    Time_RADAR = time_RADAR,
                                    FETHMM = pval.fethmm,
                                    Time_FETHMM = time_fethmm,
                                    DRME = pval.drme,
                                    Time_DRME = time_drme,
                                    TRESS = pval.TRESS,
                                    Time_TRESS = time_TRESS
                                      )
            }
          }

      
      



 if(!AddCov){
  save(AddCov, N.reps.Ctrl, N.reps.Trt, res.sim, PVAL_noCov, file = paste0("res_noCov_", sim_ind, ".rda"))
 }else{
  save(AddCov, N.reps.Ctrl, N.reps.Trt, res.sim, PVAL_Cov, file = paste0("res_Cov_", sim_ind, ".rda"))
 }



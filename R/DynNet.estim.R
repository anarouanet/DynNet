#step 1  root computation of estimation
#' Title
#'
#' @param K number of outcomes
#' @param nD number of latent processes
#' @param mapping.to.LP indicates which outcome measured which latent process, it is a mapping table between
#'  outcomes and latents processes
#' @param data indicates the data frame containing all the variables for estimating the model
#' @param if_link indicates if non linear link is used to transform an outcome
#' @param DeltaT indicates the discretization step
#' @param paras initial values for parameters
#' @param maxiter maximum iteration
#' @param nproc number of processor to be used for running this package, default value is 1
#' @param epsa threshold for the convergence criterion on the parameters, default value is 1.e-4
#' @param epsb threshold for the convergence criterion on the likelihood, default value is 1.e-4
#' @param epsd threshold for the convergence criterion on the derivatives, default value is 1.e-3
#' @param print.info to print information during the liklihood optimization, default value is FALSE
#' @param cholesky logical indicating if the variance covariance matrix is parameterized using the cholesky (TRUE, by default) or the correlation (FALSE)
#' @param MCnr number of QMC replicates to compute the integral over random effects
#' @param MCnr2 number of QMC replicates to compute the integral over random effects in the Louis Variance
#' @param nmes number of repeated measurements
#' @param data_surv dataset for survival model
#' @param predict_ui boolean indicating if bayesian estimates of random effects should be computed (FALSE by default)
#'
#' @return DynNet object
#' 
#' @import marqLevAlg randtoolbox foreach doParallel
DynNet.estim <- function(K, nD, mapping.to.LP, data, if_link = if_link, cholesky = FALSE, DeltaT=1.0, MCnr = NULL, MCnr2=NULL, nmes = NULL, data_surv = NULL, paras, 
                          maxiter = 500, nproc = 1, epsa =0.0001, epsb = 0.0001,epsd= 0.001, print.info = FALSE, predict_ui = FALSE){
  cl <- match.call()
  #  non parall Optimisation 
  # package loading
  
  debug=0

  if(debug==1 || maxiter == -1){
    
    ptm<-proc.time()

    temp <- Loglik(K = K, nD = nD, mapping =  mapping.to.LP, paraOpt = paras$paraOpt,  paraFixe = paras$paraFixe, posfix = paras$posfix, paras_k = paras$npara_k,
                   sequence = as.matrix(paras$sequence), type_int = paras$type_int, ind_seq_i = paras$ind_seq_i, MCnr = MCnr, nmes = nmes,
                   m_is = data$m_i, Mod_MatrixY = data$Mod.MatrixY, Mod_MatrixYprim = data$Mod.MatrixYprim, df=data$df,
                   x = data$x, z = data$z, q = data$q, nb_paraD = data$nb_paraD,
                   x0 = data$x0, z0 = data$z0, q0 = data$q0, cholesky = cholesky,
                   data_surv = as.matrix(data_surv), data_surv_intY = as.matrix(data$intYsurv), nYsurv = data$nYsurv, basehaz = ifelse(paras$basehaz=="Weibull", 0, 1), knots_surv = paras$knots_surv, 
                   np_surv = paras$np_surv, survival = (data$nE>0), assoc =  paras$assoc, truncation = paras$truncation, 
                   nE = data$nE, Xsurv1 = as.matrix(data$Xsurv1), Xsurv2 = as.matrix(data$Xsurv2),
                   if_link = if_link, zitr = data$zitr, ide = data$ide,
                   tau = data$tau, tau_is=data$tau_is, 
                   modA_mat = data$modA_mat, DeltaT, ii=length(data$m_i)+10)

    time=proc.time()-ptm
    h=floor(time[3]/3600)
    m=floor((time[3]-h*3600)/60)
    s=floor(time[3]-h*3600-m*60)
    cat(h,"heures ",m, "minutes ", s, "seconds","\n")
    
    para <- paras$para
    para[which(paras$posfix==0)] <- paras$paraOpt
    est <- list("fn.value" = temp)
    est$istop <- 2
    est$v <- NULL 
    est$b <- paras$paraOpt
    est$ca <- NULL 
    est$cb <- NULL
    est$rdm <- NULL
    est$ni <- -1
  }else{
    
    # marqLevAlg::marqLevAlg(b = paras$paraOpt, fn = Loglik, nproc = nproc, .packages = NULL, epsa=epsa, epsb=epsb, epsd=epsd,
    #                maxiter=maxiter, print.info = print.info, minimize = FALSE,
    #                DeltaT=DeltaT, paraFixe = paras$paraFixe, posfix = paras$posfix,
    #                paras_k = paras$npara_k, 
    #                sequence = paras$sequence, type_int = paras$type_int,
    #                K = K, nD = nD, mapping =  mapping.to.LP, m_is = data$m_i, if_link = if_link, zitr = data$zitr, ide = data$ide, 
    #                Mod_MatrixY = data$Mod.MatrixY, Mod_MatrixYprim = data$Mod.MatrixYprim, df=data$df,
    #                x = data$x, z = data$z, q = data$q, nb_paraD = data$nb_paraD,
    #                x0 = data$x0, z0 = data$z0, q0 = data$q0,tau = data$tau, tau_is=data$tau_is,
    #                modA_mat = data$modA_mat)
    # 
    # marqLevAlg::marqLevAlg(b = paras$paraOpt, fn = CInLPN:::Loglik, nproc = nproc, .packages = NULL, epsa=epsa, epsb=epsb, epsd=epsd,
    #                        maxiter=maxiter, print.info = print.info, minimize = FALSE,
    #                        DeltaT=DeltaT, paraFixe = paras$paraFixe, posfix = paras$posfix,
    #                        K = K, nD = nD, mapping =  mapping.to.LP, m_is = data$m_i, if_link = if_link,
    #                        Mod_MatrixY = data$Mod.MatrixY, Mod_MatrixYprim = data$Mod.MatrixYprim, df=data$df,
    #                        x = data$x, z = data$z, q = data$q, nb_paraD = data$nb_paraD,
    #                        x0 = data$x0, z0 = data$z0, q0 = data$q0,tau = data$tau, tau_is=data$tau_is,
    #                        modA_mat = data$modA_mat)
    # 
    #source("/Users/anais/Documents/2019 Postdoc Bordeaux/code/CInLPN/simulations_CInLPN/Thresholds/simulations/marqLevAlg_AR.R")
    #source("/Users/anais/Documents/2019 Postdoc Bordeaux/code/R/MLM/deriva_AR.R")

      ptm<-proc.time()#marqLevAlg::marqLevAlg

      temp <- try(marqLevAlg::marqLevAlg(b = paras$paraOpt, fn = Loglik, nproc = nproc, .packages = NULL, epsa=epsa, epsb=epsb, epsd=epsd,
                                         maxiter=maxiter, print.info = print.info,  minimize = FALSE,
                                         DeltaT=DeltaT, paraFixe = paras$paraFixe, posfix = paras$posfix,
                                         paras_k = paras$npara_k, 
                                         sequence = as.matrix(paras$sequence), type_int = paras$type_int, ind_seq_i = paras$ind_seq_i,  MCnr = MCnr, nmes = nmes,
                                         K = K, nD = nD, mapping =  mapping.to.LP, m_is = data$m_i, if_link = if_link, zitr = data$zitr, ide = data$ide, 
                                         Mod_MatrixY = data$Mod.MatrixY, Mod_MatrixYprim = data$Mod.MatrixYprim, df=data$df,
                                         x = data$x, z = data$z, q = data$q, nb_paraD = data$nb_paraD,
                                         x0 = data$x0, z0 = data$z0, q0 = data$q0, cholesky = cholesky, tau = data$tau, tau_is=data$tau_is,
                                         modA_mat = data$modA_mat, data_surv = as.matrix(data_surv), data_surv_intY = as.matrix(data$intYsurv), nYsurv = data$nYsurv, basehaz = ifelse(paras$basehaz=="Weibull", 0, 1), knots_surv = paras$knots_surv, 
                                         np_surv = paras$np_surv, survival = (data$nE>0), assoc =  paras$assoc, truncation = paras$truncation, 
                                         nE = data$nE, Xsurv1 = as.matrix(data$Xsurv1), Xsurv2 = as.matrix(data$Xsurv2), 
                                         clustertype="FORK", ii=length(data$m_i)+10)
                  ,silent = FALSE)
      
      time=proc.time()-ptm
      h=floor(time[3]/3600)
      m=floor((time[3]-h*3600)/60)
      s=floor(time[3]-h*3600-m*60)
      cat(h,"heures ",m, "minutes ", s, "seconds","\n")
      if(inherits(temp ,'try-error')){
        est <- list(istop=20, v=rep(0,length=((length(paras$paraOpt))*((length(paras$paraOpt)+1)/2))) ,
                    fn.value=100000000, b=paras$paraOpt, ca=1,cb=1,rdm=1,ier=-1)
      }else{
        est <- temp
      }
  }

  if(MCnr2>0){
    cat("Computation of Louis variance with ",MCnr2, " QMC draws","\n")
    N <- length(data$m_i)
    
    #library(foreach)
    #library(doSNOW)
    #library(parallel)
    
    # Si <- foreach(ii=1:2,
    #               .combine=cbind) %dopar%
    #   {
    #     test <- marqLevAlg::deriva(b = paras$paraOpt, funcpa = Loglik, nproc = 1, .packages = NULL, #epsa=epsa, epsb=epsb, epsd=epsd,
    #                                #maxiter=maxiter, print.info = print.info,  minimize = FALSE,
    #                                DeltaT=DeltaT, paraFixe = paras$paraFixe, posfix = paras$posfix,
    #                                paras_k = paras$npara_k, 
    #                                sequence = as.matrix(paras$sequence), type_int = paras$type_int, ind_seq_i = paras$ind_seq_i,  MCnr = MCnr, nmes = nmes,
    #                                K = K, nD = nD, mapping =  mapping.to.LP, m_is = data$m_i, if_link = if_link, zitr = data$zitr, ide = data$ide, 
    #                                Mod_MatrixY = data$Mod.MatrixY, Mod_MatrixYprim = data$Mod.MatrixYprim, df=data$df,
    #                                x = data$x, z = data$z, q = data$q, nb_paraD = data$nb_paraD,
    #                                x0 = data$x0, z0 = data$z0, q0 = data$q0, cholesky = cholesky, tau = data$tau, tau_is=data$tau_is,
    #                                modA_mat = data$modA_mat, data_surv = as.matrix(data_surv), data_surv_intY = as.matrix(data$intYsurv), nYsurv = data$nYsurv, basehaz = ifelse(paras$basehaz=="Weibull", 0, 1), knots_surv = paras$knots_surv, 
    #                                np_surv = paras$np_surv, survival = (data$nE>0), assoc =  paras$assoc, truncation = paras$truncation, 
    #                                nE = data$nE, Xsurv1 = as.matrix(data$Xsurv1), Xsurv2 = as.matrix(data$Xsurv2), 
    #                                #clustertype="FORK", 
    #                                ii=ii)
    #     v <- test$v[((length(paras$paraOpt)*(length(paras$paraOpt)+1))/2+1):length(test$v)]
    #     c(v)
    #   }
    
    
    if(paras$type_int == 2){
      sequence  <- sobol(n = MCnr2, dim = sum(data$q)+nD, scrambling = 1, normal = TRUE, init=T)
    }else if(paras$type_int == 1){
      sequence  <- halton(n = MCnr2, dim = sum(data$q)+nD, normal = TRUE, init=T) 
    }else if(paras$type_int == 3){
      sequence  <- torus(n = MCnr2, dim = sum(data$q)+nD, normal = TRUE, init=T) 
    }
    
    I1 <- matrix(0,length(paras$paraOpt),length(paras$paraOpt))
    I2 <- rep(0,length(paras$paraOpt))
    
    if(nproc>1){
      clustpar <- parallel::makeCluster(nproc, type="FORK")#, outfile="")
      doParallel::registerDoParallel(clustpar)    
      
      ll <- foreach(ii=1:N,
                    .combine=cbind) %dopar%
        {
          test <- marqLevAlg::deriva(b = paras$paraOpt, funcpa = Loglik, nproc = 1, .packages = NULL, #epsa=epsa, epsb=epsb, epsd=epsd,
                         #maxiter=maxiter, print.info = print.info,  minimize = FALSE,
                         DeltaT=DeltaT, paraFixe = paras$paraFixe, posfix = paras$posfix,
                         paras_k = paras$npara_k, 
                         sequence = as.matrix(sequence), type_int = paras$type_int, ind_seq_i = paras$ind_seq_i,  MCnr = MCnr2, nmes = nmes,
                         K = K, nD = nD, mapping =  mapping.to.LP, m_is = data$m_i, if_link = if_link, zitr = data$zitr, ide = data$ide, 
                         Mod_MatrixY = data$Mod.MatrixY, Mod_MatrixYprim = data$Mod.MatrixYprim, df=data$df,
                         x = data$x, z = data$z, q = data$q, nb_paraD = data$nb_paraD,
                         x0 = data$x0, z0 = data$z0, q0 = data$q0, cholesky = cholesky, tau = data$tau, tau_is=data$tau_is,
                         modA_mat = data$modA_mat, data_surv = as.matrix(data_surv), data_surv_intY = as.matrix(data$intYsurv), nYsurv = data$nYsurv, basehaz = ifelse(paras$basehaz=="Weibull", 0, 1), knots_surv = paras$knots_surv, 
                         np_surv = paras$np_surv, survival = (data$nE>0), assoc =  paras$assoc, truncation = paras$truncation, 
                         nE = data$nE, Xsurv1 = as.matrix(data$Xsurv1), Xsurv2 = as.matrix(data$Xsurv2), 
                         #clustertype="FORK", 
                         ii=ii)
          v <- test$v
          v
        }
      
      parallel::stopCluster(clustpar)
      
      for(ii in 1:N){
        I1 <- I1 + ll[,ii]%*%t(ll[,ii])
        I2 <- I2 + ll[,ii]
      }
      
    }else{
      for(ii in 1:N){
        # temp_ii <- marqLevAlg::marqLevAlg(b = paras$paraOpt, fn = Loglik, nproc = nproc, .packages = NULL, epsa=epsa, epsb=epsb, epsd=epsd,
        #                                   maxiter=maxiter, print.info = print.info,  minimize = FALSE,
        #                                   DeltaT=DeltaT, paraFixe = paras$paraFixe, posfix = paras$posfix,
        #                                   paras_k = paras$npara_k, 
        #                                   sequence = as.matrix(paras$sequence), type_int = paras$type_int, ind_seq_i = paras$ind_seq_i,  MCnr = MCnr, nmes = nmes,
        #                                   K = K, nD = nD, mapping =  mapping.to.LP, m_is = data$m_i, if_link = if_link, zitr = data$zitr, ide = data$ide, 
        #                                   Mod_MatrixY = data$Mod.MatrixY, Mod_MatrixYprim = data$Mod.MatrixYprim, df=data$df,
        #                                   x = data$x, z = data$z, q = data$q, nb_paraD = data$nb_paraD,
        #                                   x0 = data$x0, z0 = data$z0, q0 = data$q0, cholesky = cholesky, tau = data$tau, tau_is=data$tau_is,
        #                                   modA_mat = data$modA_mat, data_surv = as.matrix(data_surv), data_surv_intY = as.matrix(data$intYsurv), nYsurv = data$nYsurv, basehaz = ifelse(paras$basehaz=="Weibull", 0, 1), knots_surv = paras$knots_surv, 
        #                                   np_surv = paras$np_surv, survival = (data$nE>0), assoc =  paras$assoc, truncation = paras$truncation, 
        #                                   nE = data$nE, Xsurv1 = as.matrix(data$Xsurv1), Xsurv2 = as.matrix(data$Xsurv2), 
        #                                   clustertype="FORK", ii=ii)
        #marqLevAlg::deriva
        
        
        
        test <- deriva(b = paras$paraOpt, funcpa = Loglik, nproc = 1, .packages = NULL, #epsa=epsa, epsb=epsb, epsd=epsd,
                       #maxiter=maxiter, print.info = print.info,  minimize = FALSE,
                       DeltaT=DeltaT, paraFixe = paras$paraFixe, posfix = paras$posfix,
                       paras_k = paras$npara_k, 
                       sequence = as.matrix(sequence), type_int = paras$type_int, ind_seq_i = paras$ind_seq_i,  MCnr = MCnr2, nmes = nmes,
                       K = K, nD = nD, mapping =  mapping.to.LP, m_is = data$m_i, if_link = if_link, zitr = data$zitr, ide = data$ide, 
                       Mod_MatrixY = data$Mod.MatrixY, Mod_MatrixYprim = data$Mod.MatrixYprim, df=data$df,
                       x = data$x, z = data$z, q = data$q, nb_paraD = data$nb_paraD,
                       x0 = data$x0, z0 = data$z0, q0 = data$q0, cholesky = cholesky, tau = data$tau, tau_is=data$tau_is,
                       modA_mat = data$modA_mat, data_surv = as.matrix(data_surv), data_surv_intY = as.matrix(data$intYsurv), nYsurv = data$nYsurv, basehaz = ifelse(paras$basehaz=="Weibull", 0, 1), knots_surv = paras$knots_surv, 
                       np_surv = paras$np_surv, survival = (data$nE>0), assoc =  paras$assoc, truncation = paras$truncation, 
                       nE = data$nE, Xsurv1 = as.matrix(data$Xsurv1), Xsurv2 = as.matrix(data$Xsurv2), 
                       #clustertype="FORK", 
                       ii=ii)
        v <- test$v#[((length(paras$paraOpt)*(length(paras$paraOpt)+1))/2+1):length(test$v)]
        I1 <- I1 + v%*%t(v)
        I2 <- I2 + v
      }
    }

    V <- I1 #- 1/N*I2%*%t(I2)
    #det(-V)
    V_louis <- solve(V) #Hessian =-V, cov = [-H]^1
    #V_louis <- solve(-V)
    # (res <- Loglik(paraOpt = paras$paraOpt, DeltaT=DeltaT, paraFixe = paras$paraFixe, posfix = paras$posfix,
    #                K = K, nD = nD, mapping = mapping.to.LP, m_is = data$m_i, if_link = if_link,
    #                Mod_MatrixY = data$Mod.MatrixY, Mod_MatrixYprim = data$Mod.MatrixYprim, df=data$df,
    #                x = data$x, z = data$z, q = data$q, nb_paraD = data$nb_paraD,
    #                x0 = data$x0, z0 = data$z0, q0 = data$q0,tau = data$tau, tau_is=data$tau_is,
    #                modA_mat = data$modA_mat))
    
    #   (res <- pred(paraOpt = paras$paraOpt,paraFixe = paras$paraFixe, 
    #                  posfix = paras$posfix, DeltaT=DeltaT, K = K, nD = nD, mapping = mapping.to.LP, 
    #                  modA_mat = data$modA_mat, m_is = data$m_i, Y = data$Y, 
    #                  x = data$x, z = data$z, q = data$q, x0 = data$x0, 
    #                  z0 = data$z0, q0 = data$q0, nb_paraDu= data$nb_paraDu, 
    #                  nb_paraDw= data$nb_paraDw, tau = data$tau, tau_is=data$tau_is))    
  }


  ui=rep(0,sum(data$q)+sum(data$q0))
  #if(pred & (data$nE>0 || any(if_link==2))){
    if(predict_ui){
    ui_hat <- matrix(NA,length(data$m_i),sum(data$q)+sum(data$q0))
    maxiter=100

    for(i in 1:length(data$m_i)){
      optim_ui<- try(marqLevAlg::marqLevAlg(b = ui, paraOpt = paras$paraOpt, fn = Loglik2, nproc = nproc, .packages = NULL, epsa=epsa, epsb=epsb, epsd=epsd,
                                         maxiter=maxiter, print.info = F,  minimize = FALSE,
                                         DeltaT=DeltaT, paraFixe = paras$paraFixe, posfix = paras$posfix,
                                         paras_k = paras$npara_k, 
                                         sequence = as.matrix(paras$sequence), type_int = paras$type_int, ind_seq_i = paras$ind_seq_i,  MCnr = MCnr, nmes = nmes,
                                         K = K, nD = nD, mapping =  mapping.to.LP, m_is = data$m_i, if_link = if_link, zitr = data$zitr, ide = data$ide, 
                                         Mod_MatrixY = data$Mod.MatrixY, Mod_MatrixYprim = data$Mod.MatrixYprim, df=data$df,
                                         x = data$x, z = data$z, q = data$q, nb_paraD = data$nb_paraD,
                                         x0 = data$x0, z0 = data$z0, q0 = data$q0, cholesky = cholesky, tau = data$tau, tau_is=data$tau_is,
                                         modA_mat = data$modA_mat, data_surv = as.matrix(data_surv), data_surv_intY = as.matrix(data$intYsurv), 
                                         nYsurv = data$nYsurv, basehaz = ifelse(paras$basehaz=="Weibull", 0, 1), knots_surv = paras$knots_surv, 
                                         np_surv = paras$np_surv, survival = (data$nE>0), assoc =  paras$assoc, truncation = paras$truncation, 
                                         nE = data$nE, Xsurv1 = as.matrix(data$Xsurv1), Xsurv2 = as.matrix(data$Xsurv2), 
                                         clustertype="FORK", ii=i)
                  ,silent = T)
      if(inherits(temp ,'try-error')){
        ui_hat[i,] <- rep(NA,dim(ui_hat)[2])
      }else{
        ui_hat[i,] <- optim_ui$b
      }
    }
    est$ui_hat <- ui_hat
  }else{
    est$ui_hat <- NULL
  }
  
  
  para <- paras$para
  para[which(paras$posfix==0)] <- est$b
  
  #estimating para + fixed para
  
  est$coefficients <- para
  est$posfix <- paras$posfix
  if(MCnr2>0){
    est$LouisV <- V_louis[upper.tri(V_louis, diag=T)]
  }else{
    est$LouisV <- NULL
  }
  est
}

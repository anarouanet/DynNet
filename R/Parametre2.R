#' Function to initialize parameters k in multivariate DynNet model
#'
#' @param K number of the markers
#' @param nD number of the latent processes
#' @param nL number of exogeneous latent processes
#' @param mapping.to.LP indicates which outcome measured which latent process, it is a mapping table between
#'  outcomes and latents processes
#' @param mapping.to.LP indicates which outcome measured which exogenous latent process, it is a mapping table between
#'  outcomes and enogenous latents processes (only used for structural model)
#' @param vec_ncol_x0n vector of number of columns of model.matrix for baseline's submodel
#' @param n_col_x number of overall columns of model.matrix for change's submodel
#' @param nb_RE number of random effects
#' @param stochErr indicates if the structural model contain stochastique components
#' @param indexparaFixeUser position of parameters to be constrained
#' @param paraFixeUser values associated to the index of parameters to be constrained
#' @param L number of columns of model.matrix for temporal infuences model
#' @param paras.ini initial values for parameters, default values is NULL
#' @param ncolMod.MatrixY vector of number of columns of model.matrix for transformation submodel
#' @param link indicates link used to transform outcome
#' @param npara_k marker-specific number of transformation parameters
#' @param Survdata dataset for survival model
#' @param basehaz type of baseline hazard function
#' @param knots_surv knots for splines in baseline hazard (to develop)
#' @param assoc specification of association between longitudinal and survival models
#' @param truncation boolean for delayed entry
#' @param data dataset for longitudinal model
#' @param outcomes names of the outcomes
#' @param df degree of freedom for link==splines
#' @param nE number of survival events
#' @param np_surv cause-specific number of fixed effects in survival models
#' @param fixed.survival.models specification of survival models (without interactions)
#' @param interactionY.survival.models specification of interactions in survival models
#' @param nYsurv number of fixed effects in survival model
#' @return a list
#'
#' @importFrom stats qunif median
#' @importFrom mstate trans.comprisk msprep expand.covs
#' @importFrom survival coxph survreg
#'
#'
Parametre2 <- function(K,
                       nD,
                       nL,
                       mapping.to.LP,
                       mapping.to.LP2,
                       vec_ncol_x0n,
                       n_col_x,
                       nb_RE,
                       stochErr = FALSE,
                       indexparaFixeUser = NULL,
                       paraFixeUser = NULL,
                       L = 1,
                       paras.ini,
                       ncolMod.MatrixY,
                       link,
                       npara_k,
                       Survdata = NULL,
                       basehaz = NULL,
                       knots_surv = NULL,
                       assoc = NULL,
                       truncation = F,
                       data,
                       outcomes,
                       df,
                       nE = 0,
                       np_surv = 0,
                       fixed.survival.models = NULL,
                       interactionY.survival.models = NULL,
                       nYsurv = 0,
                       data_F = data_F) {
  cl <- match.call()
  
  #   require(MASS)
  #initialisation des parametres
  # L = number of parameters for each coefficient a of matrix A
  # K = number of outcomes
  #======================================================================================
  
  nb_paraD = nb_RE * (nb_RE + 1) / 2
  indexparaFixeForIden <- NULL
  # if user not specified initial parameters
  
  if (is.null(basehaz))
    basehaz <- "Weibull" #not to have NULL value in C++ code
  
  if (is.null(paras.ini)) {
    p <- 0 # position in the initialize parameters
    cpt1 <- 0 # compteur pour tous les paras
    cpt2 <- 0 # compteur de boucle
    #alpha_mu0
    alpha_mu0 <- rep(.1, sum(vec_ncol_x0n))
    p <- p + sum(vec_ncol_x0n)
    index_paraFixe_mu0_constraint <- NULL
    for (n in 1:nD) {
      alpha_mu0[(cpt2 + 1)] <- 0
      cpt2 <- cpt2 + vec_ncol_x0n[n]
      cpt1 <- cpt1 + vec_ncol_x0n[n]
    }
    
    #alpha_mu
    alpha_mu <- rep(.3, n_col_x)
    p <- p + n_col_x
    cpt1 <- cpt1 + n_col_x
    
    alpha_D <- rep(.1, nb_paraD)
    to_nrow <- nb_RE
    i_alpha_D <- 0
    index_paraFixeDconstraint <- NULL
    for (n in 1:nD) {
      alpha_D[i_alpha_D + 1] <- 1
      i_alpha_D <- i_alpha_D + to_nrow
      cpt1 <- cpt1 + to_nrow
      to_nrow <- to_nrow - 1
    }
    p <- p + nb_paraD
    # para of transition matrix vec_alpha_ij
    vec_alpha_ij <- rep(0.4, L * nD * nD)
    cpt1 <- cpt1 + L * nD * nD
    p <- p + L * nD * nD
    #paraB
    paraB <- NULL
    if (stochErr == TRUE) {
      paraB <- rep(.15, nD)
      cpt1 <- cpt1 + nD
      p <- p + nD
    }
    #paraSig
    paraSig <- rep(.5, K)
    cpt1 <- cpt1 + K
    p <- p + K
    ### parameters of the link function
    ParaTransformY <- c()
    for (k in 1:K) {
      if (link[k] == 'thresholds') {
        ParaTransformY <- c (ParaTransformY, c(
          2 * qunif(0.98) * (-median(data[, outcomes[k]]) + min(data[, outcomes[k]]) +
                               1) / (length(unique(data[, outcomes[k]])) - 2),
          rep(sqrt(2 * qunif(0.98) /
                     (
                       length(unique(data[, outcomes[k]])) - 2
                     )), length(unique(data[, outcomes[k]])) - 2)
        ))
        #ParaTransformY <- c(min(data[,outcomes[k]]), rep(1, length(unique(data[,outcomes[k]]))-2))
      } else{
        ParaTransformY <- c(ParaTransformY, rep(1, ncolMod.MatrixY))
      }
    }
    
    cpt1 <- cpt1 + ncolMod.MatrixY
    p <- p + ncolMod.MatrixY
    
    para_surv <- NULL
    para_basehaz <- NULL
    knots_surv <- c(0, 0) # changer !!
    #Survival
    if (!is.null(Survdata)) {
      cov_surv <- unique(unlist(strsplit(fixed.survival.models, "[+*]")))
      
      if (nE == 2) {
        tmat <- trans.comprisk(2, names = c("event-free", "event1", "event2"))
        
        Survdata$stat1 <- as.numeric(Survdata$StatusEvent == 1)
        Survdata$stat2 <- as.numeric(Survdata$StatusEvent == 2)
        Survdatalong <- msprep(
          time = c(NA, "Event", "Event"),
          status = c(NA, "stat1", "stat2"),
          data = Survdata,
          keep = cov_surv,
          trans = tmat
        )
        dim0 <- dim(Survdatalong)[2]
        Survdatalong <- expand.covs(Survdatalong, cov_surv)
        cov_survexpanded <- paste(names(Survdatalong)[(dim0 + 1):dim(Survdatalong)[2]], collapse =
                                    '+')
        form <- as.formula(
          paste(
            "Surv(time, status) ~ ",
            cov_survexpanded,
            "+ factor(trans)+ strata(trans)",
            sep = ""
          )
        )
        mod_surv <- survreg(form, data = Survdatalong, dist = "weibull")
        
        #   1/(survreg's scale)  =    rweibull shape a
        #   exp(survreg's intercept) = rweibull scale b
        # (a/b) (x/b)^(a-1) shape a, scale b
        para_basehaz <- c(
          1 / (mod_surv$coefficients[["(Intercept)"]]),
          exp(mod_surv$icoef[2]),
          1 / (mod_surv$coefficients[["(Intercept)"]] + mod_surv$coefficients[["factor(trans)2"]]),
          exp(mod_surv$icoef[3])
        )
        n1 <- 1 + np_surv[1] - 1 + ifelse(assoc %in% c(0, 1, 3, 4), 1, 2) *
          nD
        para_surv <- c(
          mod_surv$coefficients[2:(1 + np_surv[1] - 1)],
          rep(0, ifelse(assoc %in% c(0, 1, 3, 4), 1, 2) * nD),
          mod_surv$coefficients[(n1):(n1 - 1 + np_surv[1] -
                                        1)],
          rep(0, ifelse(assoc %in% c(0, 1, 3, 4), 1, 2) * nD)
        )
        
      } else if (nE == 1) {
        form <- as.formula(paste("Surv(Event, StatusEvent) ~ ", cov_surv, sep =
                                   ""))
        mod_surv <- survreg(form, data = Survdata, dist = "weibull")
        para_basehaz <- c(1 / (mod_surv$coefficients[["(Intercept)"]]), exp(mod_surv$icoef[2]))
        para_surv <- c(mod_surv$coefficients[2:(1 + np_surv - 1)], rep(0, ifelse(assoc %in%
                                                                                   c(0, 1, 3, 4), 1, 2) * nD))
      }
      p <- p + length(para_basehaz) + length(para_surv)
      np_baz <- length(para_basehaz) / nE
    }
    
    if (!is.null(interactionY.survival.models)) {
      browser()
    }
  }
  
  # if user specified initial parameters
  if (!is.null(paras.ini)) {
    map_p <- rep(NA, length(paras.ini))#attribution of parameters to latent process
    p <- 0 # position in the initialize parameters
    cpt1 <- 0 # counter for parameters
    cpt2 <- 0 # loop counter
    
    #alpha_mu0
    alpha_mu0 <- paras.ini[(p + 1):(p + sum(vec_ncol_x0n))]
    map_p[(p + 1):(p + sum(vec_ncol_x0n))] <- as.integer(sub("^LP0\\.(\\d+).*", "\\1", colnames(data_F$x0)))
    names(map_p)[(p + 1):(p + sum(vec_ncol_x0n))] <- paste0("alpha_mu0", 1:length(alpha_mu0))
    
    p <- p + sum(vec_ncol_x0n)
    index_paraFixe_mu0_constraint <- NULL
    for (n in 1:nD) {
      cpt2 <- cpt2 + vec_ncol_x0n[n]
      cpt1 <- cpt1 + vec_ncol_x0n[n]
    }
    paraFixe_mu0_constraint <- rep(1, nD)
    #alpha_mu
    alpha_mu <- paras.ini[(p + 1):(p + n_col_x)]
    map_p[(p + 1):(p + n_col_x)] <- as.integer(sub("^DeltaLP\\.(\\d+).*", "\\1", colnames(data_F$x)))
    names(map_p)[(p + 1):(p + n_col_x)] <- paste0("alpha_mu", 1:length(alpha_mu))
    p <- p + n_col_x
    cpt1 <- cpt1 + n_col_x
    #random effects
    alpha_D <- paras.ini[(p + 1):(p + nb_paraD)]
    RE_z0 <- as.integer(sub("\\(.*\\)", "", colnames(data_F$z0)))
    RE_z <- as.integer(sub("\\(.*\\)", "", colnames(data_F$z)))
    map_p[(p + 1):(p + nb_paraD)] <- rep(c(RE_z0, RE_z), times = 1:nb_RE)
    names(map_p)[(p + 1):(p + nb_paraD)] <- paste0("alpha_D", 1:length(alpha_D))
    to_nrow <- nb_RE
    i_alpha_D <- 0
    index_paraFixeDconstraint <- NULL
    
    for (n in 1:nD) {
      i_alpha_D <- i_alpha_D + to_nrow
      cpt1 <- cpt1 + to_nrow
      to_nrow <- to_nrow - 1
    }
    p <- p + nb_paraD
    
    paraFixeDconstraint <- rep(1, nD)
    # para of transition matrix vec_alpha_ij
    vec_alpha_ij <- paras.ini[(p + 1):(p + L * nD * nD)]
    map_p[(p + 1):(p + L * nD * nD)] <- rep(1:nD, each = nD)# to check if correct
    names(map_p)[(p + 1):(p + L * nD * nD)] <- paste0("vec_alpha_ij", 1:length(vec_alpha_ij))
    p <- p + L * nD * nD
    cpt1 <- cpt1 + L * nD * nD
    # paraB
    paraB <- NULL
    if (stochErr == TRUE) {
      paraB <- paras.ini[(p + 1):(p + nD)]
      map_p[(p + 1):(p + nD)] <- 1:nD
      names(map_p)[(p + 1):(p + nD)] <- paste0("paraB", 1:length(paraB))
      p <- p + nD
      cpt1 <- cpt1 + nD
    }
    #paraSig
    paraSig <- paras.ini[(p + 1):(p + K)]
    map_p[(p + 1):(p + K)] <- mapping.to.LP
    names(map_p)[(p + 1):(p + K)] <- paste0("paraSig", 1:length(paraSig))
    p <- p + K
    cpt1 <- cpt1 + K
    
    ### para of link function
    ParaTransformY <- paras.ini[(p + 1):(p + ncolMod.MatrixY)]
    i_para <- 0
    for (k in 1:K) {
      if (link[k] == "linear" & ParaTransformY[i_para + 2] == 0) {
        stop('Second parameter for linear link function cannot be set at 0 (variance)')
      }
      i_para <- i_para + npara_k[k]
    }
    
    cpt1 <- cpt1 + ncolMod.MatrixY
    map_p[(p + 1):(p + ncolMod.MatrixY)] <- rep(mapping.to.LP, times = as.numeric(table(sub(
      "\\..*", "", colnames(data_F$Mod.MatrixY)
    ))))
    names(map_p)[(p + 1):(p + ncolMod.MatrixY)] <- paste0("ParaTransformY", 1:length(ParaTransformY))
    p <- p + ncolMod.MatrixY
    
    # Weigths for formative structural model
    para_weights <- paras.ini[(p + 1):(p + nL)]
    mappingLP2LP1 <- pmin(table(mapping.to.LP2, mapping.to.LP), 1)
    mappingLP2LP1_weights <- mappingLP2LP1 * para_weights
    sum_w <- apply(mappingLP2LP1_weights, 2, sum)
    end_zero <- apply(mappingLP2LP1_weights, 2, function(x)
      any(x == 1)) & apply(mappingLP2LP1, 2, sum) > 1
    if (any(sum_w) != 1)
      stop(
        "Initial values for the weights of the formative part of the structural model need to sum to 1 within each endogenous latent process"
      )
    if (any(end_zero == TRUE))
      stop(
        "Some initial values for the weights reduce formative structure of the model by putting weight=1"
      )
    
    map_p[(p + 1):(p + nL)] <- rep(as.integer(colnames(mappingLP2LP1)), times =
                                     as.numeric(apply(mappingLP2LP1, 2, sum)))
    names(map_p)[(p + 1):(p + nL)] <- paste0("para_weights", 1:length(para_weights))
    p <- p + nL
    
    
    #Survival
    para_surv <- NULL
    para_basehaz <- NULL
    knots_surv <- c(0, 0) # changer !!
    if (!is.null(Survdata)) {
      np_baz <- ifelse(basehaz == "Weibull", 2, 0)# changer 0!!
      
      for (jj in 1:nE) {
        para_basehaz <- c(para_basehaz, paras.ini[(p + 1):(p + np_baz)])
        p <- p + np_baz  # change here?
        #}
        #for (jj in 1:nE){
        para_surv <- c(para_surv, paras.ini[(p + 1):(p + np_surv[jj])])
        p <- p + np_surv[jj] # change here?
      }
      if (basehaz == "Splines")
        cat('add number of parameters for splines in p and para_surv')
      if (basehaz == "Splines")
        cat('Define knots_surv para_basehaz')
      
    }
    
    #if(length(paras.ini) != (p + sum(df)))
    #  stop("The length of paras.ini is not correct.")
  }
  
  #final vector of initial parameters
  paras <- c(
    alpha_mu0,
    alpha_mu,
    alpha_D,
    vec_alpha_ij,
    paraB,
    paraSig,
    ParaTransformY,
    para_weights
  )
  t1 <- 0
  t2 <- 0
  
  if (nE > 0) {
    for (jj in 1:nE) {
      paras <- c(paras, para_basehaz[(t1 + 1):(t1 + np_baz)]) # change 0!!
      t1 <- t1 + np_baz
      paras <- c(paras, para_surv[(t2 + 1):(t2 + np_surv[jj])]) # change 0!!
      t2 <- t2 + np_surv[jj]
    }
  }
  
  
  if (!is.null(paras.ini)) {
    if (length(paras) != p || length(paras.ini) != p) {
      message("The length of paras.ini is not correct.")
      browser()
      stop("The length of paras.ini is not correct.")
    }
  } else{
    if (length(paras) != p)
      stop("The length of paras.ini is not correct.")
  }
  
  #initialisation
  #   paraOpt <- paras
  posfix <- rep(0, length(paras)) # 0 = non fixe 1 = fixe # initialisation
  # constraining of parameters==============
  indexFixe <- indexparaFixeForIden
  
  if (!is.null(indexparaFixeUser)) {
    if (length(indexparaFixeUser) != length(paraFixeUser)) {
      stop("The length of paraFixe does not correspond with the length of indexparaFixe")
    }
    indexFixe <- sort(unique(c(indexFixe, indexparaFixeUser)))
  }
  paraFixe <- rep(NA, length(posfix))
  if (!is.null(paraFixeUser)) {
    paraFixe[c(indexparaFixeUser)] <- paraFixeUser
  }
  paraFixe[index_paraFixe_mu0_constraint] <- rep(0, K)
  paraFixe[index_paraFixeDconstraint] <- rep(1, K)
  if (sum(!is.na(paraFixe)) == 0) {
    paraFixe = -1
    paraOpt <- paras
  } else{
    paraFixe <- paraFixe[!is.na(paraFixe)]
    posfix[indexFixe] <- 1 # fixation des paras d'indexes dans indexparaFixe
    paras[indexFixe] <- paraFixe
    paraOpt <- paras[-indexFixe]
  }
  
  if (!is.null(nL)) {
    
  }
  
  # here I start transforming parameters
  map_weights <-  map_p[which(grepl("para_weights",names(map_p)))]
  
  #alpha_mu0
  map_alpha_mu0 <- map_p[which(grepl("alpha_mu0",names(map_p)))]
  alpha_mu0_trans <- unlist(lapply(1:nD,function(x)lapply(para_weights[map_weights==x],function(l)alpha_mu0[map_alpha_mu0==x]*l)))
  # indexFixe_alpha_mu0 <- intersect(which(grepl("alpha_mu0",names(map_p))),indexFixe)
  # if(!length(indexFixe_alpha_mu)){
  #   
  # }else{
  #   indexFixe_alpha_mu_trans <- NULL
  # }
  
  #alpha_mu
  map_alpha_mu <- map_p[which(grepl("alpha_mu",names(map_p))  & !grepl("alpha_mu0",names(map_p)))]
  alpha_mu_trans <- unlist(lapply(1:nD,function(x)lapply(para_weights[map_weights==x],function(l)alpha_mu[map_alpha_mu==x]*l)))
  indexFixe_alpha_mu <- intersect(which(grepl("alpha_mu",names(map_p))  & !grepl("alpha_mu0",names(map_p))),indexFixe)
  
  #alpha_D
  # currently a nightmare so lets go on 
  # map_alpha_D <- map_p[which(grepl("alpha_D",names(map_p)))]
  # 
  # alpha_D
  # 
  # alpha_D_trans <- unlist(lapply(1:nD,function(x)lapply(para_weights[map_weights==x],function(l)alpha_D[map_alpha_D==x]*l)))
  # indexFixe_alpha_mu <- intersect(which(grepl("alpha_mu",names(map_p))  & !grepl("alpha_mu0",names(map_p))),indexFixe)
  # 
  
  #vec_alpha_ij
  vec_alpha_ij_trans <-matrix(vec_alpha_ij,nrow=nD,byrow = T)
  indexFixe_alpha_mu <- intersect(which(grepl("alpha_mu",names(map_p))  & !grepl("alpha_mu0",names(map_p))),indexFixe)
  
  
  paras_trans <- c(
    alpha_mu0_trans,
    alpha_mu_trans,
    alpha_D,
    vec_alpha_ij,
    paraB,
    paraSig,
    ParaTransformY
  )
  
  
  
  
  
  return(
    list(
      para = paras_trans,
      paraOpt = paraOpt,
      paraFixe = paraFixe,
      posfix = posfix,
      L = L,
      basehaz = basehaz,
      knots_surv = knots_surv,
      np_surv = np_surv,
      assoc = assoc,
      truncation = truncation
    )
  )
}

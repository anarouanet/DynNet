#' function that returns the parameters object required in the DynNet function.
#' 
#' This function creates the vectors of initial parameter values (paras.ini), 
#' the indices of user-specified fixed parameters (Fixed.para.index), 
#' and their corresponding fixed values (Fixed.para.values), which are required to 
#' run the DynNet function.
#' 
#' @param structural.model a list of 7 arguments used to specify the structural model: \describe{
#'    \item{\code{fixed.LP0}}{a one-sided linear formula object for specifying the 
#'    fixed effects in the submodel for the baseline level of processes. Note that 
#'    there is no need to specify a random effect model for the baseline level of 
#'    processes as we systematically set a process-specific random intercept (with 
#'    variance fixed to 1 for identifiability purpose). For identifiability purposes, 
#'    the mean intercepts are fixed to 0 (not estimated).}
#' 
#'    \item{\code{fixed.DeltaLP}}{a two-sided linear formula object for 
#'    specifying the response outcomes (one the left part of ~ symbol) 
#'    and the covariates with fixed-effects (on the right part of ~ symbol) 
#'    in the submodel for change over time of latent processes.}
#' 
#'    \item{\code{random.DeltaLP}}{a one-sided linear formula object 
#'    for specifying the random effects in the submodel for change 
#'    over time of latent processes.}
#' 
#'    \item{\code{trans.matrix}}{a one-sided linear formula object 
#'    for specifying a model for elements of the transition matrix, 
#'    which captures the temporal influences between latent processes.}
#' 
#'    \item{\code{fixed.survival}}{a one-sided linear formula object 
#'    for specifying the covariates in the survival sub-model. In competing 
#'    risks model, the specification for the two events should be separated 
#'    by the "|" symbol.}
#' 
#'    \item{\code{interactionY.survival}}{a one-sided linear formula 
#'    object for specifying the covariates in interaction with the 
#'    dynamics of the latent processes, in the survival sub-model. 
#'    In competing risks model, the specification for the two events 
#'    should be separated by the "|" symbol. Only additional terms 
#'    should be included (No "*" symbol). Covariates in 
#'    \code{interactionY.survival} should also be included in \code{fixed.survival.}}
#' 
#'    \item{\code{delta.time}}{indicates the discretisation step to be used for latent processes}
#'    }
#' 
#' @param measurement.model is a list of arguments detailed below used to specify the measurement model: \describe{
#' 
#' \item{\code{link}}{ indicates the link functions to be used to transform the outcomes. It takes values in "linear" for a linear transformation and "n-type-d" for a I-splines transformation where "n" indicates the number of nodes, "type" (which takes values in "quant", "manual", "equi") indicates where the nodes are placed, and "d" indicates the degree of the I-splines.}
#' 
#' \item{\code{knots}}{ argument indicates if necessary the place of knots (when placed manually with "manual"), default value is NULL}
#' }
#' 
#' @param Time indicates the name of the covariate representing the time 
#' @param subject indicates the name of the covariate representing the grouping structure
#' @param data indicates the data frame containing all the variables for estimating the model.
#' @param Tentry name of the variable of entry time
#' @param Event name of the variable of event time
#' @param StatusEvent name of the variable of event status
#' @param basehaz type of baseline hazard function
#' @param assocT specifies the type of association between the time-to-event(s) and the latent process(es). Values include "r.intercept" for random intercept, "r.slope" for random slope, "r.intercept/slope" for both random intercept and slope, "c.value" for current value. 
#' @param p.initlev initial values for the fixed effects of the model describing the initial level of the processes
#' @param p.slope initial values for the fixed effects of the model describing the slope of the processes.
#' @param varcovRE initial values for the lower-triangular matrix in the Cholesky decomposition of the random-effects variance–covariance matrix.
#' @param transitionmatrix initial values for the transition matrix
#' @param var.errors initial values for the marker-specific measurement error variances.
#' @param transformationY initial values for the marker-specific transformation functions (see links in DynNet help)
#' @param baseline1 initial values for the baseline hazard function in the time-to-event model (in the survival setting) or in the first transition model (in the competing risks setting).
#' @param p.X1 initial values for the covariate fixed effects (excluding association and latent-process interaction parameters)in the time-to-event model (in the survival setting) or in the first transition model (in the competing risks setting).
#' @param fix.p.X1 indicator if the parameters \code{p.X1} are fixed (e.g. not estimated).
#' @param p.asso.int1 initial values for the interactions between covariates and functions of the latent processes in the time-to-event model (in the survival setting) or in the first transition model (in the competing risks setting).
#' @param fix.p.asso.int1 indicator if the parameters \code{p.asso.int1} are fixed (e.g. not estimated).
#' @param baseline2 initial values for the baseline hazard function in the second transition model (competing risks setting only). Default to NULL.
#' @param fix.baseline2 indicator if the parameters \code{baseline2} are fixed (e.g. not estimated).
#' @param p.X2 initial values for the covariate fixed effects (excluding association and latent-process interaction parameters) in the second transition model (competing risks setting only). Default to NULL.
#' @param fix.p.X2 indicator if the parameters \code{p.X2} are fixed (e.g. not estimated).
#' @param p.asso2 initial values for the association parameters in the second transition model (competing risks setting only). Default to NULL.
#' @param fix.p.asso2 indicator if the parameters \code{p.asso2} are fixed (e.g. not estimated).
#' @param p.asso.int2 initial values for the interactions between covariates and functions of the latent processes in the second transition model (competing risks setting only). Default to NULL.
#' @param fix.p.asso.int2 indicator if the parameters \code{p.asso.int2} are fixed (e.g. not estimated).
#' @return A list with the following elements:
#' \describe{
#'   \item{paras.ini}{Vector of initial values for all parameters.}
#'   \item{Fixed.para.index}{Vector containing the indexes of the parameters to be fixed.}
#'   \item{Fixed.para.values}{Vector of fixed values corresponding to the parameters specified in \code{Fixed.para.index}.}
#'}
#' @export
#' @importFrom splines bs
#' @examples
#' ### example 1
#' Delta <- 1
#' structural.model = list(fixed.LP0 = ~ 1 + C2 | 1 + C2,
#'                                       fixed.DeltaLP = L2 | L3  ~ 1 | 1 ,
#'                                       random.DeltaLP = ~ 1|1,
#'                                       trans.matrix = ~ 1,
#'                                       delta.time = Delta)
#' measurement.model = list(link.functions = list(links = c(NULL,NULL),
#'                                       knots = list(NULL, NULL)))
#'                                       
#' Parameters <- enter_param(structural.model = structural.model,
#'               measurement.model = measurement.model,
#'               Time = "time",
#'               subject = "id",
#'               data = data,
#'               p.initlev = c(0.000, 0.059, 0.000, 0.163), 
#'               fix.p.initlev = c(1, 0, 1, 0), 
#'               p.slope = c(-0.050, -0.153), 
#'               varcovRE =  c(1.000, 0.000, 0.322, 0.000, 1.000, 0.000, 0.077, 0.139, 0.000, 0.177), 
#'               fix.varcovRE = c(1, 1, 0, 1, 1, 1, 0, 0, 1, 0),
#'               transitionmatrix = c(-0.354, 0.114, 0.116, -0.090), 
#'               var.errors = c(0.287, 0.554), 
#'               transformationY = c(1.107, 1.889, 0.881, 1.329) )
#'
#' mod1 <- DynNet(structural.model = structural.model,
#'               measurement.model = measurement.model,
#'               parameters = Parameters,
#'               option = list(nproc = 1, print.info = TRUE,  MCnr = 10, 
#'                             univarmaxiter = 7, epsa = 1e-5, epsb = 1e-4, epsd = 1e-2),
#'               Time = "time",
#'               subject = "id",
#'               data = data
#' )
#' 
#' summary(mod1)
#' 
#'  ### example 2
#'  Delta <- 0.5
#'  structural.model = list(fixed.LP0 = ~ 1 + C2 | 1 + C2,
#'                          fixed.DeltaLP = L1 + L2| L3  ~ 1 + time| 1 + time,
#'                          random.DeltaLP = ~ 1|1,
#'                          trans.matrix = ~ 1 + bs(x = time, knots =c(2), 
#'                                                  intercept = F, degree = 2),
#'                          delta.time = Delta)
#' measurement.model = list(link.functions = list(links = c(NULL,NULL, NULL),
#'                                                    knots = list(NULL, NULL, NULL)))
#'  Parameters <- enter_param(structural.model = structural.model,
#'                            measurement.model = measurement.model,
#'                            Time = "time",
#'                            subject = "id",
#'                            data = data,
#'                            p.initlev = c(0.000,0.065, 0.000, 0.168), 
#'                            fix.p.initlev = c(1, 0, 1, 0), 
#'                            p.slope = c(-0.054, 0.000, -0.119, -0.009), 
#'                            varcovRE =  c(1.000, 0.000, 0.473, 0.000, 1.000, 0.000, 0.057, -0.182,
#'                                          0.000, 0.174), 
#'                            fix.varcovRE = c(1, 1, 0, 1, 1, 1, 0, 0, 1, 0),
#'                            transitionmatrix = c(-0.523, 0.000, 0.000, 0.000, 0.073, 0.083, 
#'                                                 0.119, 0.106, 0.015, 0.157, 0.054, 0.087, 
#'                                                 -0.079, 0.000, 0.000, 0.000), 
#'                            fix.transitionmatrix = c(0,rep(1,3),rep(0,9),rep(1,3)),
#'                            var.errors = c(8.876, 0.287, 0.546), 
#'                            transformationY = c(0.531, 0.256, 1.100,
#'                                                1.891,  0.846, 1.345))
#'                                                
#' mod2 <- DynNet(structural.model = list(fixed.LP0 = ~ 1 + C2 | 1 + C2,
#'                                      fixed.DeltaLP = L1 + L2| L3  ~ 1 + time| 1 + time,
#'                                      random.DeltaLP = ~ 1|1,
#'                                      trans.matrix = ~ 1 + bs(x = time, knots =c(2), 
#'                                                              intercept = F, degree = 2),
#'                                      delta.time = Delta),
#'              measurement.model = list(link.functions = list(links = c(NULL,NULL, NULL),
#'                                                             knots = list(NULL, NULL, NULL))),
#'              
#'              parameters = Parameters,
#'              option = list(nproc = 2, print.info = TRUE, MCnr = 10, 
#'                            univarmaxiter = 7, epsa = 1e-5, epsb = 1e-4, epsd = 1e-2),
#'              Time = "time",
#'              subject = "id",
#'              data = data
#' )
#'
#' summary(mod2)
#'
#'
#'       
enter_param<-function(structural.model, 
                      measurement.model,
                      Time,
                      subject,
                      data,
                      Tentry = NULL,
                      Event = NULL,
                      StatusEvent = NULL,
                      basehaz = NULL,
                      assocT=NULL,
                      p.initlev,
                      fix.p.initlev=rep(0,length(p.initlev)),
                      p.slope,
                      fix.p.slope=rep(0,length(p.slope)),
                      varcovRE,
                      fix.varcovRE=rep(0,length(varcovRE)),
                      transitionmatrix,
                      fix.transitionmatrix=rep(0,length(transitionmatrix)),
                      parab=NULL,
                      fix.parab=rep(0,length(parab)),
                      var.errors,
                      fix.var.errors=rep(0,length(var.errors)),
                      transformationY,
                      fix.transformationY=rep(0,length(transformationY)),
                      baseline1=NULL,
                      fix.baseline1=rep(0,length(baseline1)),
                      p.X1=NULL,
                      fix.p.X1=rep(0,length(p.X1)),
                      p.asso1=NULL,
                      fix.p.asso1=rep(0,length(p.asso1)),
                      p.asso.int1=NULL,
                      fix.p.asso.int1=rep(0,length(p.asso.int1)),
                      baseline2=NULL,
                      fix.baseline2=rep(0,length(baseline2)),
                      p.X2=NULL,
                      fix.p.X2=rep(0,length(p.X2)),
                      p.asso2=NULL,
                      fix.p.asso2=rep(0,length(p.asso2)),
                      p.asso.int2=NULL,
                      fix.p.asso.int2=rep(0,length(p.asso.int2))
){
  ### check if all component of the model specification are well filled ####
  if(missing(structural.model))stop("The argument structural.model must be specified")
  if(missing(measurement.model))stop("The argument measurement.model must be specified")
  #if(missing(parameters))stop("The argument parameters must be specified")
  
  if(is.null(structural.model$fixed.DeltaLP))stop("The argument structural.model$fixed.DeltaLP must be specified")
  if(is.null(structural.model$trans.matrix))stop("The argument structural.model$trans.matrix must be specified")
  
  
  if(is.null(measurement.model$link.functions) || all(is.null(measurement.model$link.functions$links))){
    links <- NULL
    knots <- NULL
    measurement.model$link.functions =list(links = links, knots = knots)
  }else{
    if(all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="quant", x))) ==0)&& 
       all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="manual", x))) ==0)&& 
       all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="equi", x))) ==0)&& 
       all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="linear", x))) ==0)&& 
       all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="thresholds", x))) ==0))
      stop("The only available link functions are 'linear', 'splines' and 'thresholds' functions.")
  }
  
  #if(is.null(parameters$Fixed.para.index))stop("The argument parameters$Fixed.para.index cannot be NULL")
  #if(is.null(parameters$Fixed.para.values))stop("The argument parameters$Fixed.para.values cannot be NULL")
  
  
  survival= FALSE
  if(!is.null(structural.model$fixed.survival)){
    survival = TRUE
    fixed.survival <- structural.model$fixed.survival
  }
  
  ### identification of model components #####
  ## components of structural model
  fixed_X0 <- structural.model$fixed.LP0
  fixed_DeltaX <- structural.model$fixed.DeltaLP
  randoms_DeltaX <- structural.model$random.DeltaLP
  mod_trans <- structural.model$trans.matrix
  
  interactionY.survival <- NULL
  if(!is.null(structural.model$interactionY.survival)){
    interactionY.survival <- structural.model$interactionY.survival
  }
  
  # components of measurement model
  link <- measurement.model$link.functions$links
  knots <- measurement.model$link.functions$knots
  ## components of parameters initialisation
  #indexparaFixeUser <- parameters$Fixed.para.index
  
  #paraFixeUser <- parameters$Fixed.para.values
  #paras.ini <- parameters$paras.ini
  ## component of option
  
  #### fixed effects pre-traitement ####
  
  ### for DeltaLP
  # if(missing(fixed_DeltaX)) stop("The argument fixed_DeltaX must be specified in any model")
  if(!inherits(fixed_DeltaX,"formula")) stop("The argument fixed_DeltaX must be a formula")
  
  ### outcomes and latent processes ####
  outcome <- as.character(attr(terms(fixed_DeltaX),"variables"))[2]
  outcomes_by_LP<-strsplit(outcome,"[|]")[[1]]
  nD <- length(outcomes_by_LP) # nD: number of latent process
  
  outcomes <- NULL
  mapping.to.LP <- NULL
  for(n in 1:nD){
    outcomes_n <- strsplit(outcomes_by_LP[n],"[+]")[[1]]
    outcomes_n <-as.character(sapply(outcomes_n,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE))
    outcomes_n <- unique(outcomes_n)
    if(is.null(outcomes_n)) stop("at least one marker must be specified for each latent process" )
    outcomes <- c(outcomes, outcomes_n)
    mapping.to.LP <- c(mapping.to.LP, rep(n,length(outcomes_n)))
  }
  
  K <- length(outcomes)
  all.Y<-seq(1,K)
  
  fixed_DeltaX.model=strsplit(gsub("[[:space:]]","",as.character(fixed_DeltaX)),"~")[[3]]
  fixed_DeltaX.models<-strsplit(fixed_DeltaX.model,"[|]")[[1]]# chaque model d'effet fixe mais en vu de connaitre tous les pred.fixed du modele multi
  
  if(nD !=length(fixed_DeltaX.models)) stop("The number of models does not correspond to the number of latent processes")
  
  if(nD > K){
    stop("There are too many latent processes compared to the indicated number of markers")
  }
  
  ### pre-traitement of fixed effect on initial levels of processes
  if(is.null(fixed_X0)){
    fixed_X0<- ~1
    fixed_X0.models <- rep("1",nD)
  }
  if(!inherits(fixed_X0,"formula")) stop("The argument fixed_X0 must be a formula")
  
  fixed_X0.models =strsplit(gsub("[[:space:]]","",as.character(fixed_X0)),"~")[[2]]
  fixed_X0.models<- as.vector(strsplit(fixed_X0.models,"[|]")[[1]]) 
  for(nd in 1:nD){
    if(fixed_X0.models[nd]=="~-1") fixed_X0.models[nd] <-"~1" # au moins l'intcpt
  }
  
  ### pre-traitement of fixed effect on survival 
  colnames<-colnames(data)
  fixed.survival.models <- NULL
  assoc <- 0
  if(survival){
    if(is.null(fixed.survival)){
      fixed.survival<- ~1
    }else{
      if(!inherits(fixed.survival,"formula")) stop("The argument fixed.survival must be a formula") 
    }
    
    fixed.survival.models <- strsplit(gsub("[[:space:]]","",as.character(fixed.survival)),"~")[[2]]
    covsurv <- unique(as.vector(strsplit(fixed.survival.models,"[|*+]")[[1]]))
    fixed.survival.models <- as.vector(strsplit(fixed.survival.models,"[|]")[[1]]) 
    if(!all(covsurv[-which(covsurv=="1")]%in%colnames))stop("All covariates in fixed.survival should be in the dataset")
    
    if(!is.null(assocT)){
      if(!assocT%in%c("r.intercept", "r.slope", "r.intercept/slope", "c.value"))
        stop("assocT should be defined as r.intercept, r.slope, r.intercept/slope, c.value.")
      assoc <- switch(assocT, "r.intercept"=0, "r.slope"=1, "r.intercept/slope"=2, "c.value"=3, "c.slope"=4, "c.value/slope"=5)
    }else{
      assocT <- "r.intercept/slope"
      assoc <- 2 # random intercept and slope
    }
    
    if(assoc <= 2){
      message(" add interactions ui * X in survival model")
      
    }
  }
  
  ### pre-traitement of interactions with Y on survival 
  
  interactionY.survival.models <- NULL
  if(survival){
    
    if(!is.null(interactionY.survival)){
      if(!inherits(interactionY.survival,"formula")) stop("The argument interactionY.survival must be a formula") 
      interactionY.survival.models <- strsplit(gsub("[[:space:]]","",as.character(interactionY.survival)),"~")[[2]]
      intYsurv <- (as.vector(strsplit(interactionY.survival.models,"[|*+]")[[1]]))
      interactionY.survival.models <- as.vector(strsplit(interactionY.survival.models,"[|]")[[1]]) 
      if(!all(intYsurv%in%colnames))stop("All covariates in interactionY.survival should be in the dataset")
      if(any(grepl( "*", interactionY.survival, fixed = TRUE)))stop("Only + terms should be included in interactionY.survival, no *.")
    }
  }
  
  
  ### pre-traitement of random effect on processes  intercept and slope
  #### randoms effet on DeltaLP 
  randoms_X0.models <- rep("1",nD)
  #### randoms effet on DeltaX
  
  if(missing(randoms_DeltaX) || is.null(randoms_DeltaX)){
    randoms_DeltaX<- ~1
    randoms_DeltaX.models <- rep("1",nD)
  }
  if(!inherits(randoms_DeltaX,"formula")) stop("The argument random must be a formula")
  randoms_DeltaX.model=strsplit(gsub("[[:space:]]","",as.character(randoms_DeltaX)),"~")[[2]]
  randoms_DeltaX.models<-strsplit(randoms_DeltaX.model,"[|]")[[1]]    
  
  #### traitement of  mod_trans: transition matrix##### 
  if(missing(mod_trans)){
    mod_trans <- ~ 1 # constant transition matrix
  } 
  if(!inherits(mod_trans,"formula")) stop("The argument mod_trans must be a formula")
  mod_trans.model=strsplit(gsub("[[:space:]]","",as.character(mod_trans)),"~")[[2]]
  
  if(nD!=length(fixed_X0.models)){
    stop("The number of models for initial latent processes does not correspond with the number of latent processes")
  }
  if(nD!=length(fixed_DeltaX.models)){
    stop("The number of models for the change over time of latent processes does not correspond with the number of latent processes")
  }
  
  ### traitement of transformation models ##
  if(is.null(link)){
    link <- rep("linear",K)
  }else if(length(link)!=K) stop("The number transformation links must be equal to the number of markers")
  
  if(any(link == "thresholds")){
    j     <- which(link == 'thresholds')
    
    Y0    <- data[,outcomes[j]]
    minY0 <- apply(as.matrix(Y0), 2, min, na.rm=TRUE)# min(Y0,rm.na=TRUE)
    maxY0 <- apply(as.matrix(Y0), 2, max, na.rm=TRUE)# max(Y0,rm.na=TRUE)
    
    if(length(j==1))
      ide0 <- matrix(0, nrow = length(j),  ncol = max(sapply(j, function(x) max(Y0, na.rm =T)-min(Y0, na.rm =T) ))) #change dimensions
    else
      ide0 <- matrix(0, nrow = length(j),  ncol = max(sapply(j, function(x) max(Y0[,x], na.rm =T)-min(Y0[,x], na.rm =T) ))) #change dimensions
    
    nbzitr0 <- 2
    #idlink0 <- 3
    #ntrtot0 <- as.integer(maxY0 - minY0)
    if (!(all(is.integer(minY0) | !all(is.integer(maxY0)))))
      stop("With the thresholds link function, the longitudinal outcome must be discrete")
    
    zitr <- c()
    for(i in 1:length(j)){
      if(length(j)>1){
        Y0tmp <- Y0[,i]
      }else{
        Y0tmp <- Y0
      }
      
      if (!all(Y0tmp[which(!is.na(Y0tmp))] %in% minY0[i]:maxY0[i]))
        stop("With the threshold link function, problem with the outcome data, must be discrete")
      
      IND <- sort(unique(Y0tmp))
      IND <- IND[1:(length(IND) - 1)] - minY0[i] + 1
      #ide0 <- rep(0, as.integer(maxY0[i] - minY0[i])) #change dimensions
      ide0[i, IND] <- 1
      
      #if(i==1)
      #  zitr <- matrix(rep(0, length(j)*nbzitr0), length(j), nbzitr0)
      #zitr[i, nbzitr0] <- maxY0[i]
      zitr <- c(zitr, minY0[i], maxY0[i])
    }
    
    if(!is.null(type_int)){
      if(!type_int %in% c("MC", "sobol", "halton", "torus"))
        stop("With the thresholds link function, type_int should be either antithetic, sobol, halton or torus. antithetic not developed yet, sorry.")
    }
    
  }else{
    zitr <- 0
    ide0  <- 0 
  }
  
  
  ## If joint model -  Event and StatusEvent data
  if(survival){
    
    basehaz <- ifelse(!is.null(basehaz), basehaz, "Weibull")
    
    if(!basehaz%in%c("Weibull",'splines'))
      stop('basehaz should be either Weibull or splines.')
  }
  
  
  #DataFormat==============================================================
  

  id_and_Time <- data[,c(subject,Time)]
  Ni<-unique(data[,subject])
  I <- length(Ni) # number of visit
  K <- length(outcomes)
  
  # Pre-traitement of data : delete lignes with no observation
  d <- as.data.frame(data[,outcomes])
  R <- as.numeric(apply(X = d, MARGIN = 1, FUN = is_na_vec))
  data <-data[which(R==0),]
  
  #all predictor of the model==============================================================
  all.pred.fixed_X0 <- NULL
  all.pred.fixed_DeltaX <- NULL
  all.pred.randoms_X0 <- NULL
  all.pred.randoms_DeltaX <- NULL
  all.pred.mod_trans <- NULL
  for( n in 1: nD){
    all.pred.fixed_X0 <- c(all.pred.fixed_X0,list(strsplit(fixed_X0.models[n],"[+]")[[1]]))
    all.pred.fixed_DeltaX <- c(all.pred.fixed_DeltaX,list(strsplit(fixed_DeltaX.models[n],"[+]")[[1]]))
    all.pred.randoms_X0 <- c(all.pred.randoms_X0,list(strsplit(randoms_X0.models[n],"[+]")[[1]]))
    all.pred.randoms_DeltaX <- c(all.pred.randoms_DeltaX,list(strsplit(randoms_DeltaX.models[n],"[+]")[[1]]))
  }
  #
  all.pred.fixed_X0<-sapply(all.pred.fixed_X0,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
  all.pred.fixed_X0<-sapply(all.pred.fixed_X0,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
  #
  all.pred.fixed_DeltaX<-sapply(all.pred.fixed_DeltaX,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
  all.pred.fixed_DeltaX<-sapply(all.pred.fixed_DeltaX,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
  #
  all.pred.randoms_X0<-sapply(all.pred.randoms_X0,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
  all.pred.randoms_X0<-sapply(all.pred.randoms_X0,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
  #
  all.pred.randoms_DeltaX<-sapply(all.pred.randoms_DeltaX,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
  all.pred.randoms_DeltaX<-sapply(all.pred.randoms_DeltaX,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
  #
  all.pred.mod_trans <- c(all.pred.mod_trans,list(strsplit(mod_trans.model,"[+]")[[1]]))
  all.pred.mod_trans<-sapply(all.pred.mod_trans,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
  all.pred.mod_trans<-sapply(all.pred.mod_trans,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
  #ajout
  all.preds<-unlist(unique(c(unlist(all.pred.fixed_X0), unlist(all.pred.fixed_DeltaX), 
                             unlist(all.pred.randoms_X0), unlist(all.pred.randoms_DeltaX),
                             all.pred.mod_trans)))
  all.preds<-inclu_intercerpt(all.preds) # remplace 1 par "(Intercept)"
  all.preds<-unique(unlist(sapply(all.preds,FUN=function(x)strsplit(x,":")[[1]])))
  
  ###all.pred_san_inter signifie all.pred_sans_intercept
  all.preds<-all.preds[-which(all.preds %in% c("(Intercept)"))]
  all.preds <- c(all.preds,Time)
  all.preds <- unique(all.preds[which(all.preds %in% colnames)])  
  
  #Case of  unobserved components  at time t
  m_i <-as.data.frame(table(as.factor(data[,subject])))$Freq # matrice of frequencies m_i
  DeltaT <- structural.model$delta.time
  tau_is <- data[,Time]/DeltaT # vector of individuals visits vectors
  tau_is <- as.numeric(as.character(tau_is)) # verify that all integer?
  
  
  Tmax <- max(tau_is,na.rm = TRUE)
  
  Survdata <- NULL
  if(survival){
    if(!(Event%in%colnames)) stop("Event should be in the dataset")
    if(!(StatusEvent%in%colnames)) stop("StatusEvent should be in the dataset")
    
    basehaz <- ifelse(!is.null(basehaz), basehaz, "Weibull")
    if(Tentry != "Tentry" & !(Tentry%in%colnames)) stop("Tentry should be in the dataset")
    
    if(!basehaz%in%c("Weibull",'splines'))
      stop('basehaz should be either Weibull or splines.')
    
    if(basehaz == "Splines")
      cat("Define knots_surv here")
    #One survival dataframe with one line per individual
    first_line <- sapply(unique(data[,subject]), function(x) which(data[,subject]==x)[1])
    
    if(!(Tentry %in%names(data))) data$Tentry <- 0
    Survdata <- data[first_line, c(Tentry, Event, StatusEvent)]
    names(Survdata) <- c("Tentry", "Event", "StatusEvent")
    
    if(length(which(covsurv!="1"))>0){
      Survdata <- cbind(Survdata, data[first_line, covsurv])
      names(Survdata)[(dim(Survdata)[2]-length(covsurv)+1):dim(Survdata)[2]]<-covsurv
    }
  }
  if(survival && assoc %in%c(3,5))
    Tmax <- max(Tmax, round(max(Survdata$Event,na.rm = TRUE)/DeltaT))#If
  
  if(!survival)
    Survdata<-NULL
  
  tau <- 0:Tmax  
  Y <- NULL
  IND <- NULL
  indY <- NULL
  data0 <- NULL
  
  ###creation de data0==========
  ## for x and z
  all.Y<-seq(1,K)
  for (k in 1:K)
  {
    dtemp <- data[,c(subject,outcomes[k],all.preds)]
    Y <- c(Y, dtemp[,outcomes[k]])
    IND <- c(IND, dtemp[,subject])
    indY <- c(indY,rep(all.Y[k],nrow(dtemp)))
    data0 <- rbind(data0, dtemp[,c(setdiff(colnames(dtemp),outcomes[k]))])
  }
  
  data0<-cbind(data0,Y,indY)
  data0<-data0[order(data0[,subject]),]
  data0<-na.omit(data0)
  ## for x and z
  indY <- data0[,"indY"]
  Y<-data0[,"Y"]
  
  
  #### only for x0 and x ###############
  xs <- as.data.frame(unique(data0[,c(subject,setdiff(all.preds,Time))])) # this does not work with time-dependent covariates
  colnames(xs) <- c(subject,setdiff(all.preds,Time))
  qsz <- as.data.frame(xs[rep(row.names(xs), rep(length(tau), dim(xs)[1])),])
  colnames(qsz) <- c(subject,setdiff(all.preds,Time))
  Times <- as.data.frame(DeltaT*rep(tau, I))
  colnames(Times) <- Time
  
  if(dim(qsz)[1] != dim(Times)[1]){
    stop("Covariates (other than time) should be time-independant.")
  }
  
  data_c0 <- cbind(qsz,Times)
  data_xzMatA_cov <-data_c0  
  data_xzMatA_cov <-data_xzMatA_cov[order(data_xzMatA_cov[,subject], data_xzMatA_cov[,Time]),]
  data_xzMatA_cov <- data_xzMatA_cov[!duplicated(data_xzMatA_cov),]
  x_cov<- NULL
  for(n in 1:nD){
    indLP_x <- rep(n, dim(data_xzMatA_cov)[1])
    data_x_cov_i <- cbind(data_xzMatA_cov, indLP_x)
    x_cov <- rbind(x_cov, data_x_cov_i)
  }
  x_cov <- x_cov[order(x_cov[,subject],x_cov[,Time]),]
  
  ##only for x0 #####
  x0 <- NULL
  nb_x0_n <- NULL
  col_n<-list()
  x0_cov <- x_cov[which(x_cov[,Time]==0),]
  #
  ##
  indLP_x0 <- x0_cov$indLP_x
  for(n in 1:nD){
    r<-as.formula(paste(subject, fixed_X0.models[n], sep="~-1+"))
    x0n<-model.matrix(r,data=x0_cov)
    nb_x0_n <- c(nb_x0_n,ncol(x0n))
    if(length(x0n)==0){
      col <- paste(n,"zero",sep="")
      x0n<-matrix(assign(col,rep(0,dim(x0_cov)[1])))
      nb_x0_n <- c(nb_x0_n,ncol(x0n))
    }
    
    colnames<-colnames(x0n)
    colnames<-paste("LP0",n,colnames,sep=".")
    colnames(x0n) <-colnames
    col_n <-c(col_n,list(colnames))
    x0<-cbind(x0,as.matrix(x0n))
  }
  x0 <-cbind(indLP_x0,x0)
  ### remplissage avec les zeros
  tous_col_x0 <-unlist(col_n)
  for(i in 1:nrow(x0)){
    col_i <- unlist(col_n[[x0[i,"indLP_x0"]]])
    col_0<-tous_col_x0[which(!(tous_col_x0 %in% col_i))]
    x0[i,col_0]<-0 #  passer pour optimisation
  }
  x0 <- as.matrix(x0)
  colnames <- colnames(x0)
  x0 <- as.matrix(x0[,-c(1)])
  colnames(x0) <- colnames[-c(1)]
  #   x0 <- as.matrix(x0)
  
  ##only for x #####
  x <- NULL
  nb_x_n <- NULL
  col_n<-list()
  indLP_x <- x_cov$indLP_x
  for(n in 1:nD){
    r<-as.formula(paste(subject, fixed_DeltaX.models[n], sep="~-1+"))
    xn<-model.matrix(r,data=x_cov)
    nb_x_n <- c(nb_x_n,ncol(xn))
    if(length(xn)==0){
      col <- paste(n,"zero",sep="")
      xn<-matrix(assign(col,rep(0,dim(x_cov)[1])))
      nb_x_n <- c(nb_x_n,ncol(xn))
    }
    colnames<-colnames(xn)
    colnames<-paste("DeltaLP",n,colnames,sep=".")
    colnames(xn) <-colnames
    col_n <-c(col_n,list(colnames))
    x<-cbind(x,as.matrix(xn))
  }
  x <-cbind(indLP_x,x)
  ### filling with zeros
  tous_col_x <-unlist(col_n)
  for(i in 1:nrow(x)){
    col_i <- unlist(col_n[[x[i,"indLP_x"]]])
    col_0<-tous_col_x[which(!(tous_col_x %in% col_i))]
    x[i,col_0]<-0 # z  passer pour optimisation
  }
  
  x <- as.matrix(x)
  colnames <- colnames(x)
  x <- as.matrix(x[,-c(1)])
  colnames(x) <- colnames[-c(1)]
  
  #===================================================================
  #     construction  of matrices z0 et z========================
  data_z_cov <- data_xzMatA_cov[, c(subject,Time)]
  z_cov<- NULL
  for(n in 1:nD){
    indY_z <- rep(n, dim(data_z_cov)[1])
    data_z_cov_i <- cbind(data_z_cov, indY_z)
    z_cov <- rbind(z_cov, data_z_cov_i)
  }
  z_cov <- z_cov[order(z_cov[,subject],z_cov[,Time]),]
  #   z_cov[,Time] <- z_cov[,Time]*DeltaT ######
  
  #### only for z0 ####
  z0_cov <- z_cov[which(z_cov[,Time]==0),]
  indY_z0 <- z0_cov$indY_z
  z0 <- NULL
  col_n<-list()
  q0 <- NULL
  nb_paraDw <- 0
  for(n in 1:nD){
    r<-as.formula(paste(subject,randoms_X0.models[n], sep="~-1+"))
    z0n<-model.matrix(r,data=z0_cov)
    if(length(z0n)==0){
      col <- paste(n,"zero",sep="")
      z0n<-matrix(assign(col,rep(0,dim(z0_cov)[1])))
    }
    colnames<-colnames(z0n)
    colnames<-paste(n,colnames,sep="")
    colnames(z0n) <-colnames
    col_n <-c(col_n,list(colnames))
    z0<-cbind(z0,z0n)
    q0 <- c(q0,ncol(z0n))
  }
  
  z0 <-cbind(indY_z0,z0)
  ### filling with zeros
  tous_col_z0 <-unlist(col_n)
  for(i in 1:nrow(z0)){
    col_i <- unlist(col_n[[z0[i,"indY_z0"]]])
    col_0<-tous_col_z0[which(!(tous_col_z0 %in% col_i))]
    z0[i,col_0]<-0 # z  passer pour optimisation
  }
  z0 <- z0[,-c(1)]
  z0 <- as.matrix(z0)
  
  
  #### only for z ####
  indY_z <- z_cov$indY_z
  z <- NULL
  col_n<-list()
  q <- NULL
  nb_paraDu <- 0
  
  for(n in 1:nD){
    r<-as.formula(paste(subject,randoms_DeltaX.models[n], sep="~-1+"))
    zn<-model.matrix(r,data=z_cov)
    if(length(zn)==0){
      col <- paste(n,"zero",sep="")
      zn<-matrix(assign(col,rep(0,dim(z_cov)[1])))
    }
    colnames<-colnames(zn)
    colnames<-paste(n,colnames,sep="")
    colnames(zn) <-colnames
    col_n <-c(col_n,list(colnames))
    z<-cbind(z,zn)
    q <- c(q,ncol(zn))
  }
  
  if(all(randoms_DeltaX.models=="-1"))
    q <- rep(0,nD)
  if(length(which(randoms_DeltaX.models=="-1"))>0 & length(which(randoms_DeltaX.models=="-1"))<length(randoms_DeltaX.models))
    stop('If one dimension does not have a random slope, the other dimensions should not either.')
  
  z <-cbind(indY_z,z)
  ### filling with zeros
  tous_col_z <-unlist(col_n)
  for(i in 1:nrow(z)){
    col_i <- unlist(col_n[[z[i,"indY_z"]]])
    col_0<-tous_col_z[which(!(tous_col_z %in% col_i))]
    z[i,col_0]<-0 #
  }
  z <- z[,-c(1)]
  z <- as.matrix(z)
  #============================================================
  # design matrix for transition model
  f<-as.formula(paste(subject,mod_trans.model, sep="~"))# put subject, just to have a left side for the formula
  modA_mat<-model.matrix(as.formula(paste(subject,mod_trans.model, sep="~")),data=data_xzMatA_cov)
  
  #   dim(modA_mat)
  #   head(modA_mat)
  #============================================================
  #design matrix for markers transformation
  Y <- as.matrix(data[,outcomes])
  tr_Y <- f.link(outcomes = outcomes, Y=as.data.frame(Y), link=link, knots =knots)
  Mod.MatrixY <- tr_Y$Mod.MatrixY
  Mod.MatrixYprim <- tr_Y$Mod.MatrixYprim
  knots <- tr_Y$knots
  degree = tr_Y$degree
  df <- tr_Y$df
  minY <- tr_Y$minY
  maxY <- tr_Y$maxY
  
  nb_RE <- sum(q0,q)
  nb_paraD <- nb_RE*(nb_RE+1)/2
  
  #If joint model
  Event <- NULL
  StatusEvent <- NULL
  
  if(!is.null(Survdata)){ # Interesting for development to multi-state and interval censoring...
    type=ifelse(length(unique(Survdata[,3]))>2, "mstate", "right")
    surv_obj <- survival::Surv(Survdata[,2], Survdata[,3], type=type)
    Event <- surv_obj[,1]
    StatusEvent <- surv_obj[,2]
    #message("check use of mstate here....")
  }
  
  nE <- length(fixed.survival.models)
  Xsurv1 <- 0
  Xsurv2 <- 0
  np_surv <- NULL
  intYsurv <- 0
  nYS <- rep(0,2)
  
  if(nE>0){
    for(n in 1: nE){
      all.pred.fixed.survival.models <- list(strsplit(fixed.survival.models[n],"[+*]")[[1]])
      Xsurv <- as.matrix(model.matrix(as.formula(paste("",fixed.survival.models[n], sep="~")),data=Survdata)[,-1])
      
      if(n==1)
        Xsurv1 <- Xsurv
      if(n==2)
        Xsurv2 <- Xsurv
      np_surv <- c(np_surv, dim(Xsurv)[2] + ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD)
    }   
    nYS <- rep(0,nE) #to have a vector in cpp program
    intYsurv <- NULL
    
    if(!is.null(interactionY.survival.models)){
      
      tmp_intYS <- as.vector(strsplit(interactionY.survival.models,"[|]")[[1]])
      for(j in 1:nE){
        intYsurv<- cbind(intYsurv, as.matrix(model.matrix(as.formula(paste("",interactionY.survival.models[j], sep="~")),data=Survdata)[,-1]))
        
        tmp_intYS_j<- as.vector(strsplit(tmp_intYS[j],"[+]")[[1]])
        nYS[j]<-length(tmp_intYS_j)
      }
      np_surv <- np_surv + nYS*ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD
    }else{
      intYsurv <- 0
    }
  }else{
    np_surv <-0
  }
  
  vec_ncol_x0n <- nb_x0_n # number of parameters on initial level of processes
  n_col_x <- ncol(x) # number of parameters on processes slope
  L <- ncol(modA_mat)
  ncolMod.MatrixY <- ncol(Mod.MatrixY)
  
  
  assocT <- NULL
  if(!is.null(assoc)){
    assocT <- ifelse(assoc==0, "r.intercept",ifelse(assoc==1, "r.slope",ifelse(assoc==2, "r.intercept/slope",ifelse(
      assoc==3, "c.value",ifelse(assoc==4, "c.slope","c.value/slope")
    ))))
  }
  
  ### Start Parametre.R ###
  npara_k <- sapply(outcomes, function(x) length(grep(x, names(data.frame(Mod.MatrixY)))))

  nb_paraD = nb_RE*(nb_RE+1)/2
  indexparaFixeForIden <- NULL
  # if user not specified initial parameters
  
  if(is.null(basehaz))
    basehaz<-"Weibull" #not to have NULL value in C++ code
  
  # if user specified initial parameters
  
  p <- 0 # position in the initialize parameters
  cpt1 <-0 # counter for parameterd
  cpt2<-0 # loop counter
  #alpha_mu0

  if(length(p.initlev)!=length((p+1):(p+sum(vec_ncol_x0n)))){
    stop("p.initlev should contain ",length((p+1):(p+sum(vec_ncol_x0n)))," parameters.")
  }
  if(length(fix.p.initlev)!=length((p+1):(p+sum(vec_ncol_x0n)))){
    stop("fix.p.initlev should contain ",length((p+1):(p+sum(vec_ncol_x0n)))," parameters.")
  }

  
  alpha_mu0 <- p.initlev#paras.ini[(p+1):(p+sum(vec_ncol_x0n))]
  p <- p+ sum(vec_ncol_x0n)
  index_paraFixe_mu0_constraint <-NULL
  for(n in 1:nD){
    #alpha_mu0[(cpt2+1)] <- 0
    cpt2 <- cpt2 + vec_ncol_x0n[n]
    cpt1 <- cpt1 + vec_ncol_x0n[n]
  }
  paraFixe_mu0_constraint <- rep(1,nD)
  #alpha_mu
  if(length(p.slope)!=length((p+1):(p+n_col_x))){
    stop("p.slope should contain ",length((p+1):(p+n_col_x))," parameters.")
  }
  if(length(fix.p.slope)!=length((p+1):(p+n_col_x))){
    stop("fix.p.slope should contain ",length((p+1):(p+n_col_x))," parameters.")
  }
  alpha_mu <- p.slope#paras.ini[(p+1):(p+n_col_x)]
  p <- p+n_col_x
  cpt1 <- cpt1 + n_col_x
  
  #alpha_D parameters for cholesky of all random effects
  if(length(varcovRE)!=length((p+1):(p+nb_paraD))){
    stop("varcovRE should contain ",length((p+1):(p+nb_paraD))," parameters.")
  }
  if(length(fix.varcovRE)!=length((p+1):(p+nb_paraD))){
    stop("fix.varcovRE should contain ",length((p+1):(p+nb_paraD))," parameters.")
  }
  
  alpha_D <- varcovRE#paras.ini[(p+1):(p+nb_paraD)]
  to_nrow <- nb_RE
  i_alpha_D <- 0
  index_paraFixeDconstraint <- NULL
  
  for(n in 1:nD){
    #if(link[n] != "thresholds")
    #alpha_D[i_alpha_D+1] <- 1
    i_alpha_D <- i_alpha_D + to_nrow
    cpt1 <- cpt1 + to_nrow
    to_nrow <- to_nrow -1
  }
  p <- p+nb_paraD
  paraFixeDconstraint <- rep(1,nD)
  # para of transition matrix vec_alpha_ij
  #alpha_D parameters for cholesky of all random effects

  if(length(transitionmatrix)!=length((p+1):(p + L*nD*nD))){
    stop("transitionmatrix should contain ",length((p+1):(p + L*nD*nD))," parameters.")
  }
  if(length(fix.transitionmatrix)!=length((p+1):(p + L*nD*nD))){
    stop("fix.transitionmatrix should contain ",length((p+1):(p + L*nD*nD))," parameters.")
  }
  vec_alpha_ij <- transitionmatrix#paras.ini[(p+1):(p + L*nD*nD)]
  p <- p + L*nD*nD
  cpt1 <- cpt1 + L*nD*nD
  # paraB
  paraB <- NULL
  stochErr = F
  if(stochErr==TRUE){
    if(length(parab)!=length((p+1):(p + nD))){
      stop("paraB should contain ",length((p+1):(p + nD))," parameters.")
    }
    if(length(fix.parab)!=length((p+1):(p + nD))){
      stop("fix.paraB should contain ",length((p+1):(p + nD))," parameters.")
    }
    paraB <- parab#paras.ini[(p+1):(p + nD)]
    p <- p + nD
    cpt1 <- cpt1 + nD
  }
  #paraSig
  if(length(var.errors)!=length((p+1):(p + K))){
    stop("var.errors should contain ",length((p+1):(p + K))," parameters.")
  }
  if(length(fix.var.errors)!=length((p+1):(p + K))){
    stop("fix.var.errors should contain ",length((p+1):(p + K))," parameters.")
  }
  paraSig <- var.errors#paras.ini[(p+1):(p + K)]
  p <- p + K
  cpt1 <- cpt1 + K
  
  ### para of link function
  if(length(transformationY)!=length((p+1):(p + ncolMod.MatrixY))){
    stop("transformationY should contain ",length((p+1):(p + ncolMod.MatrixY))," parameters.")
  }
  if(length(fix.transformationY)!=length((p+1):(p + ncolMod.MatrixY))){
    stop("fix.transformationY should contain ",length((p+1):(p + ncolMod.MatrixY))," parameters.")
  }
  ParaTransformY <- transformationY#paras.ini[(p+1):(p + ncolMod.MatrixY)]
  i_para <- 0
  for(k in 1:K){
    if(link[k]=="linear" & ParaTransformY[i_para+2]==0){
      stop('Second parameter for linear link function cannot be set at 0 (variance)')
    }
    i_para <- i_para + npara_k[k]
  }
  
  cpt1 <- cpt1 + ncolMod.MatrixY
  p <- p + ncolMod.MatrixY
  
  
  #Survival
  para_surv <- NULL
  para_basehaz <- NULL
  knots_surv <- c(0,0) # changer !!
  if(!is.null(Survdata)){
    # if(nE ==1){
    #   np_surv <- dim(Survdata)[2]-3 + ifelse(assoc%in%c(0, 1, 3, 4),1,2)
    # }else{
    #   np_surv <- dim(Survdata)[2]-3 + ifelse(assoc%in%c(0, 1, 3, 4),1,2)
    # }
    np_baz <- ifelse(basehaz=="Weibull",2, 0)# changer 0!!

    for (jj in 1:nE){
      if(jj==1){
        if(length(baseline1)!=length((p+1) : (p + np_baz)))
          stop("baseline1 should contain ",length((p+1) : (p + np_baz))," parameters.")
        if(length(fix.baseline1)!=length((p+1) : (p + np_baz)))
          stop("fix.baseline1 should contain ",length((p+1) : (p + np_baz))," parameters.")
        para_basehaz <- c(para_basehaz, baseline1)#paras.ini[(p+1) : (p + np_baz)])  
        
      }else{
        
        if(length(baseline2)!=length((p+1) : (p + np_baz)))
          stop("baseline2 should contain ",length((p+1) : (p + np_baz))," parameters.")
        if(length(fix.baseline2)!=length((p+1) : (p + np_baz)))
          stop("fix.baseline2 should contain ",length((p+1) : (p + np_baz))," parameters.")
        para_basehaz <- c(para_basehaz, baseline2)#paras.ini[(p+1) : (p + np_baz)])  
        
      }
      #para_basehaz <- c(para_basehaz, paras.ini[(p+1) : (p + np_baz)])  
      p <- p + np_baz  # change here?
      #}
      #for (jj in 1:nE){
      if(jj==1){
        
        if(length(p.X1)!=dim(Xsurv)[2])
          stop("p.X1 should contain ",dim(Xsurv)[2]," parameters.")
        if(length(fix.p.X1)!=dim(Xsurv)[2])
          stop("fix.p.X1 should contain ",dim(Xsurv)[2]," parameters.")
        
        if(length(p.asso1)!=ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD)
          stop("p.asso1 should contain ",ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD," parameters.")
        if(length(fix.p.asso1)!=ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD)
          stop("fix.p.asso1 should contain ",ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD," parameters.")
        
        if(length(p.asso.int1)!=nYS[1]*ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD)
          stop("p.asso.int1 should contain ",nYS[1]*ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD," parameters.")
        if(length(fix.p.asso.int1)!=nYS[1]*ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD)
          stop("fix.p.asso.int1 should contain ",nYS[1]*ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD," parameters.")
        
        
        
        
        para_surv <- c(para_surv, p.X1, p.asso1, p.asso.int1)#paras.ini[(p+1) : (p + np_baz)])  
        
        
      }else{
        
        if(length(p.X2)!=dim(Xsurv)[2])
          stop("p.X2 should contain ",dim(Xsurv)[2]," parameters.")
        if(length(fix.p.X2)!=dim(Xsurv)[2])
          stop("fix.p.X2 should contain ",dim(Xsurv)[2]," parameters.")
        
        if(length(p.asso2)!=ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD)
          stop("p.asso2 should contain ",ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD," parameters.")
        if(length(fix.p.asso2)!=ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD)
          stop("fix.p.asso2 should contain ",ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD," parameters.")
        
        if(length(p.asso.int2)!=nYS[2]*ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD)
          stop("p.asso.int2 should contain ",nYS[2]*ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD," parameters.")
        if(length(fix.p.asso.int2)!=nYS[2]*ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD)
          stop("fix.p.asso.int2 should contain ",nYS[2]*ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD," parameters.")

        para_surv <- c(para_surv, p.X2, p.asso2, p.asso.int2)#paras.ini[(p+1) : (p + np_baz)])  
        
      }
      #para_surv <- c(para_surv, paras.ini[(p + 1 ) : (p + np_surv[jj])]) 
      p <- p + np_surv[jj] # change here?
    }
    if(basehaz=="Splines") cat('add number of parameters for splines in p and para_surv')
    if(basehaz=="Splines") cat('Define knots_surv para_basehaz')
    
  }
  #if(length(paras.ini) != (p + sum(df)))
  #  stop("The length of paras.ini is not correct.")
  
  
  #final vector of initial parameters
  paras <- c(alpha_mu0, alpha_mu, alpha_D, vec_alpha_ij,  paraB, paraSig, ParaTransformY)
  t1 <- 0
  t2 <- 0
  
  if(nE>0){
    for(jj in 1:nE){
      paras <- c(paras, para_basehaz[(t1+1) : (t1 + np_baz)]) # change 0!!
      t1 <- t1 + np_baz
      paras <- c(paras, para_surv[(t2 + 1) : (t2 + np_surv[jj])]) # change 0!!
      t2 <- t2 + np_surv[jj]
    }
  }

  # if(nE>0){
  #   for(jj in 1:nE){
  #     paras <- c(paras, para_basehaz[(t1+1) : (t1 + np_baz)]) # change 0!!
  #     t1 <- t1 + np_baz
  #   }
  # 
  #   for(jj in 1:nE){
  #     paras <- c(paras, para_surv[(t2 + 1) : (t2 + np_surv[jj])]) # change 0!!
  #     t2 <- t2 + np_surv[jj]
  #   }
  # }
  indparaFixeUser <- c(fix.p.initlev, fix.p.slope, fix.varcovRE, fix.transitionmatrix, 
                         fix.parab, fix.var.errors, fix.transformationY, fix.baseline1, 
                         fix.p.X1, fix.p.asso1, fix.p.asso.int1, fix.baseline2, fix.p.X2, 
                         fix.p.asso2, fix.p.asso.int2)
  if(!all(indparaFixeUser==0)){
    indexparaFixeUser <- which(indparaFixeUser==1)
    Fixed.para.values <- paras[indexparaFixeUser]
  }else{
    stop("at least one parameter should be fixed for identifiability purposes.")
  }
  
  return(list("paras.ini"=paras, "Fixed.para.index"=indexparaFixeUser, "Fixed.para.values"=Fixed.para.values))
}
  

  
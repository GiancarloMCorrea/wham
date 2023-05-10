set_LAA = function(input, LAA, growth)
{
  data = input$data
  par = input$par
  map = input$map
    
  n_re_par = 4 # number of parameters for RE

  # LAA default options:
  LAA_re_ini = matrix(0, ncol = data$n_ages, nrow = data$n_years_model)
  data$LAA_re_model = 1 # none default
  data$LAA_est = rep(0, data$n_ages)
  data$SDLAA_est = c(0,0)
  data$isG_nonparametric = 0L # default: it is not nonparametric
  LAA_ini = log( 60 + (8 - 60)*exp(-0.2*(1:data$n_ages - 1)) )
  SD_ini_LAA = c(log(3), log(7)) # CV1 and CVA

  if(!is.null(LAA)) {

    data$isG_nonparametric = 1L # it is nonparametric
    if(is.null(growth)) data$isG_parametric = 0L 

    if(!is.null(LAA$LAA_vals)) LAA_ini = log(LAA$LAA_vals)
    if(!is.null(LAA$est_pars)) data$LAA_est[LAA$est_pars] = 1
    if(!is.null(LAA$re))  {
      if(!(LAA$re %in% c("none","iid","iid_a","ar1_a","2dar1","3dar1"))) stop("LAA$re must be one of the following: 'none','iid','iid_a','ar1_a','2dar1','3dar1'")
      data$LAA_re_model <- match(LAA$re, c("none","iid","iid_a","ar1_a","2dar1","3dar1")) # Respect this order to create array later
    }

    if(!is.null(LAA$SD_vals)){
      if(length(LAA$SD_vals) != 2) stop("Number of 'SD_vals' must be equal to 2.")
      SD_ini_LAA <- log(LAA$SD_vals)
    }

    if(!is.null(LAA$SD_est)){
        if(length(LAA$SD_est) > 2) stop("Number of 'SD_est' should be equal or less than 2.")
        data$SDLAA_est[LAA$SD_est] = 1
    }

  }

  # This section only used for 3D smoother:
  data$ay3D_IndexG = as.matrix(expand.grid("age" = seq_len(data$n_ages), "year" = seq_len(data$n_years_model) ))
  if(is.null(data$Var3D_ParamG)) data$Var3D_ParamG = 0 # Var_Param == 0 Conditional, == 1 Marginal


  # growth pars --------------------------
  
  par$LAA_a = LAA_ini
  par$LAA_re = LAA_re_ini
  par$LAA_repars = rep(0, times = n_re_par)
  par$LAA_repars[1] = log(0.1) # start sigma at 0.1, rho at 0
  par$SDLAA_par = SD_ini_LAA

  # --------------------------------
  # Map LAA:
  # Default map LAA
  map$LAA_a <- factor(rep(NA, times = data$n_ages))
  map$SDLAA_par <- factor(rep(NA, times = 2))

  if(is.null(growth)) { # only do it if nonparametric
    tmp1 <- par$LAA_a
    tmp1[data$LAA_est==0] = NA
    ind.notNA <- which(!is.na(tmp1))
    tmp1[ind.notNA] <- 1:length(ind.notNA)
    map$LAA_a <- factor(tmp1)

    # SD pars:
    tmp1 <- par$SDLAA_par
    tmp1[data$SDLAA_est==0] = NA
    ind.notNA <- which(!is.na(tmp1))
    tmp1[ind.notNA] <- 1:length(ind.notNA)
    map$SDLAA_par <- factor(tmp1)
  }

  # RE info:
  # LAA_re: "none","iid","iid_a","ar1_a","2dar1"
  tmp <- par$LAA_re
  if(data$LAA_re_model == 1) tmp[] = NA # no RE (either estimate RE for all ages or none at all)
  if(data$LAA_re_model %in% c(2,5,6)){ # iid 2dar1 and 3dar1
    tmp[] = 1:(dim(tmp)[1]*dim(tmp)[2]) # all y,a estimated
  }
  if(data$LAA_re_model %in% c(3,4)){ # ar1_a (devs by age, constant by year)
    for(i in 1:dim(tmp)[2]) tmp[,i] = i
  }
  map$LAA_re <- factor(tmp)

  # LAA_repars: 
  if(data$LAA_re_model == 1) tmp <- rep(NA,n_re_par) # no RE pars to estimate
  if(data$LAA_re_model == 2) tmp <- c(1,NA,NA,NA) # estimate sigma
  if(data$LAA_re_model == 3) tmp <- c(1,NA,NA,NA) # estimate sigma over ages
  if(data$LAA_re_model == 4) tmp <- c(1,2,NA,NA) # ar1_a: estimate sigma, rho_a
  if(data$LAA_re_model == 5) tmp <- c(1:3,NA) # 2dar1: estimate all
  if(data$LAA_re_model == 6) tmp <- c(1:4) # 3dar1: estimate all
  map$LAA_repars = factor(tmp)

  # End section
  input$data = data
  input$par = par
  input$map = map
  return(input)

}
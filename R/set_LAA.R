set_LAA = function(input, LAA)
{
  data = input$data
  par = input$par
  map = input$map
  
  n_RE_cols = 4 # For RE parameters
  # sigma, rho_a, rho_y, rho_c

  # LAA default options:
  n_par_def = 3 
  # 3 parameters (vonB): K, Linf, L1
  # 4 parameters (Richards): K, Linf, L1, gamma
  # n_ages parameters (nonparametric): 1, ..., n_ages
  data$LAA_model = 1 # 1: vB-classic, 2: Richards, 3: Nonparametric
  data$LAA_est <- rep(0, times = n_par_def) # default = don't estimate LAA parameters
  data$SD_est = c(0,0) # SD parameters estimate?
  data$n_LAA_par = n_par_def
  data$n_LAA_dim = 3
  # 3 RE dim (vonB): K, Linf, L1
  # 4 RE dim (Richards): K, Linf, L1, gamma
  # 1 RE dim (nonparametric)
  data$LAA_re_model = rep(1, times = data$n_LAA_dim) # default = no RE / 'none'
  LAA_re_ini = array(0, dim = c(data$n_years_model, data$n_ages, data$n_LAA_dim))
  LAA_ini = c(log(0.2), log(60), log(8)) # K, Linf, and L1
  SD_ini = c(log(3), log(7)) # SD1 and SDA
  n_par = c(3, 4, data$n_ages) # number of parameters for each method (vB: 3, and Richards: 4, nonparametric: n_ages)
  n_RE_dim = c(3, 4, 1) # number of RE model for each method
  
  # prepare LAA options:
  if(!is.null(LAA)){
  
	# Some mandatory inputs should be provided:
	if(is.null(LAA$model)) stop("LAA$model must be provided.")
	if(is.null(LAA$init_vals)) stop("LAA$init_vals must be provided.")

	# LAA model:
    #if(!is.null(LAA$model)){ # LAA model to be used
    if(!(LAA$model %in% c("vB_classic", "Richards", "Nonparametric"))) stop("LAA$model must be 'vB_classic', 'Richards', or 'Nonparametric'")
    data$LAA_model <- match(LAA$model, c("vB_classic", "Richards", "Nonparametric")) 
    data$n_LAA_par = n_par[data$LAA_model] # number of parameters to estimate
    data$LAA_est = rep(0, times = n_par[data$LAA_model]) # estimate?
    data$n_LAA_dim = n_RE_dim[data$LAA_model]
	data$LAA_re_model = rep(1, times = data$n_LAA_dim) # default = no RE / 'none'
    LAA_re_ini = array(0, dim = c(data$n_years_model, data$n_ages, n_RE_dim[data$LAA_model]))
    #}

	# LAA re model:
    if(!is.null(LAA$re)){
      if(length(LAA$re) != data$n_LAA_dim) stop("Length of 'LAA$re' must be equal to the number of RE dimensions of the chosen LAA model.")
      if(data$LAA_model %in% c(1,2)) {
		  for(k in 1:data$n_LAA_par) {
			  if(!(LAA$re[k] %in% c("none","iid_y","iid_c","ar1_y","ar1_c"))) stop(paste0("LAA$re[", k, "] must be one of the following: 'none','iid_y','iid_c','ar1_y','ar1_c'"))
			  data$LAA_re_model[k] <- match(LAA$re[k], c("none","iid_y","iid_c","ar1_y","ar1_c")) # Respect this order to create array later
		  }
	  }
	  if(data$LAA_model %in% c(3)) {
		if(!(LAA$re %in% c("none","iid","2dar1","3dgmrf"))) stop("LAA$re must be one of the following: 'none','iid','2dar1','3dgmrf'")
		data$LAA_re_model <- match(LAA$re, c("none","iid","2dar1","3dgmrf")) # Respect this order to create array later			
	  }
    }
	
	# Init vals: mandatory
    #if(!is.null(LAA$init_vals)){
    if(length(LAA$init_vals) != data$n_LAA_par) stop("Length of 'LAA$init_vals' must be equal to the number of parameters of the chosen LAA model.")
    LAA_ini <- log(LAA$init_vals)
    #}

    if(!is.null(LAA$SD_vals)){
      if(length(LAA$SD_vals) != 2) stop("Length of 'LAA$SD_vals' must be 2.")
      SD_ini <- log(LAA$SD_vals)
    }
  
    if(!is.null(LAA$est_pars)){
        if(length(LAA$est_pars) > data$n_LAA_par) stop("Length of 'LAA$est_pars' should be equal or less than the number of parameters of the chosen LAA model.")
        data$LAA_est[LAA$est_pars] = 1
    }

    if(!is.null(LAA$SD_est)){
        if(length(LAA$SD_est) > 2) stop("Length of 'LAA$SD_est' should be equal or less than 2.")
        data$SD_est[LAA$SD_est] = 1
    }

  } 
  
  #data$n_LAA_est <- sum(data$LAA_est)
  # This section only used for 3D smoother (Nonparametric):
  data$ay3D_IndexL = as.matrix(expand.grid("age" = seq_len(data$n_ages), "year" = seq_len(data$n_years_model) ))
  if(is.null(data$Var3D_ParamL)) data$Var3D_ParamL = 0 # Var_Param == 0 Conditional, == 1 Marginal


  # Update pars --------------------------
  
  par$LAA_a = LAA_ini
  par$LAA_re = LAA_re_ini
  par$LAA_repars = matrix(0, ncol = n_RE_cols, nrow = data$n_LAA_dim)
  par$LAA_repars[,1] = log(0.1) # start sigma at 0.1, rho at 0
  par$SD_par = SD_ini

  # --------------------------------
  # Map LAA:
  
  # Main pars:
  tmp1 <- par$LAA_a
  tmp1[data$LAA_est==0] = NA
  ind.notNA <- which(!is.na(tmp1))
  tmp1[ind.notNA] <- 1:length(ind.notNA)
  map$LAA_a <- factor(tmp1)  

  # SD pars:
  tmp1 <- par$SD_par
  tmp1[data$SD_est==0] = NA
  ind.notNA <- which(!is.na(tmp1))
  tmp1[ind.notNA] <- 1:length(ind.notNA)
  map$SD_par <- factor(tmp1)

  # RE info:
  map$LAA_re = NULL
  map$LAA_repars = NULL
  tmp.g.repars <- matrix(NA, ncol = n_RE_cols, nrow = data$n_LAA_dim)
  
  # Parametric equations:
  if(data$LAA_model %in% c(1,2)) {
    max_val_par = 0
	active_sum = -1
	this_max = 0
	for(i in 1:data$n_LAA_par) {
		tmp1 <- matrix(NA, nrow = data$n_years_model, ncol = data$n_ages)
		if(data$LAA_re_model[i] %in% c(2,4)){ # iid ar1- only y
		  tmp1[] = rep(1:nrow(tmp1), times = ncol(tmp1))  # all y estimated
		  active_sum = active_sum + 1
		  max_val_par = max_val_par + this_max * min(1, active_sum)
		  this_max = max(tmp1, na.rm = TRUE)
		}
		if(data$LAA_re_model[i] %in% c(3,5)){ # iid ar1 - only c
		  loop_row = rep(0, times = ncol(tmp1) + nrow(tmp1) - 1)
		  loop_row[1:ncol(tmp1)] = (ncol(tmp1) - 1):0
		  loop_col = rep(ncol(tmp1) - 1, times = ncol(tmp1) + nrow(tmp1) - 1)
		  loop_col[(length(loop_col) - ncol(tmp1) + 1):length(loop_col)] = (ncol(tmp1) - 1):0

		  for(j in seq_along(loop_col)) {
			tmp1[(loop_row[j]:loop_col[j])*(nrow(tmp1) + 1) + (j - ncol(tmp1) + 1)] <- j
		  }
		  active_sum = active_sum + 1
		  max_val_par = max_val_par + this_max * min(1, active_sum)
		  this_max = max(tmp1, na.rm = TRUE)
		}

		map$LAA_re <- c(map$LAA_re, as.vector(tmp1) + max_val_par)
		# LAA_repars: 
		if(data$LAA_re_model[i] == 1) tmp.g.repars[i,] <- rep(NA,n_RE_cols) # no RE pars to estimate
		if(data$LAA_re_model[i] == 2) tmp.g.repars[i,] <- c(1,NA,NA,NA) # estimate sigma y
		if(data$LAA_re_model[i] == 3) tmp.g.repars[i,] <- c(1,NA,NA,NA) # estimate sigma c
		if(data$LAA_re_model[i] == 4) tmp.g.repars[i,] <- c(1,NA,2,NA) # ar1_y: estimate sigma y, rho_y
		if(data$LAA_re_model[i] == 5) tmp.g.repars[i,] <- c(1,NA,2,NA) # ar1_c: estimate sigma c, rho_c
	}
	map$LAA_re = as.factor(map$LAA_re)
	ind.notNA <- which(!is.na(tmp.g.repars))
	tmp.g.repars[ind.notNA] <- 1:length(ind.notNA)
	map$LAA_repars = factor(tmp.g.repars)
  }
  
  # Nonparametric:
  if(data$LAA_model %in% c(3)) { 
	  tmp <- par$LAA_re
	  if(data$LAA_re_model == 1) tmp[] = NA # no RE (either estimate RE for all ages or none at all)
	  if(data$LAA_re_model %in% 2:4){ # iid 2dar1 and 3dar1
		tmp[] = 1:(dim(tmp)[1]*dim(tmp)[2]) # all y,a estimated
	  }
	  map$LAA_re <- factor(tmp)
	  # LAA_repars: 
	  if(data$LAA_re_model == 1) tmp.g.repars <- rep(NA,n_RE_cols) # no RE pars to estimate
	  if(data$LAA_re_model == 2) tmp.g.repars <- c(1,NA,NA,NA) # estimate sigma
	  if(data$LAA_re_model == 3) tmp.g.repars <- c(1:3,NA) # 2dar1: estimate all
	  if(data$LAA_re_model == 4) tmp.g.repars <- c(1:4) # 3dar1: estimate all
	  map$LAA_repars = factor(tmp.g.repars)
  }
  
  # End section
  input$data = data
  input$par = par
  input$map = map
  return(input)

}
set_WAA = function(input, waa_opts = NULL, WAA)
{
  # Empirical Weight-at-age data ----------------------------------------------------------------------------------

  asap3 = if(is.null(input$asap3)) asap3 = NULL
  else asap3 = input$asap3
  
  data = input$data
  data$use_catch_waa = matrix(0, nrow = data$n_years_model, ncol = data$n_fleets)
  data$use_index_waa = matrix(0, nrow = data$n_years_model, ncol = data$n_indices)
  data$isW_parametric = 0L # default value
  if(!is.null(asap3)) {
	  i <- c(seq(1,(asap3$n_fleets+1)*2-1,2),(asap3$n_fleets+1)*2 + 1:2)
	  WAA_pointers <- asap3$WAA_pointers[i] #wham has no discard data, so remove those WAA matrices
	  data$waa = array(NA, dim = c(length(asap3$WAA_mats), data$n_years_model, data$n_ages))
	  for(i in 1:length(asap3$WAA_mats)) data$waa[i,,] = asap3$WAA_mats[[i]]
	  data$waa_pointer_indices = asap3$index_WAA_pointers
      data$waa_pointer_fleets = WAA_pointers[1:data$n_fleets]
      data$waa_pointer_totcatch = WAA_pointers[data$n_fleets + 1]
      data$waa_pointer_ssb = WAA_pointers[data$n_fleets + 2]
      data$waa_pointer_jan1 = WAA_pointers[data$n_fleets + 3]
      # if ASAP3 input provided, automatically use EWAA
	  data$waa_cv = array(0, dim = dim(data$waa))
	} else {
		if(is.null(waa_opts[["waa"]])){
			if(is.null(waa_opts$waa_pointer_totcatch) | is.null(waa_opts$waa_pointer_ssb) | is.null(waa_opts$waa_pointer_jan1) | is.null(waa_opts$waa_pointer_fleets) | is.null(waa_opts$waa_pointer_indices)) {
			  stop("If 'waa' is not available, waa pointers must be provided: 
				'waa_pointer_totcatch', 'waa_pointer_ssb', 'waa_pointer_jan1', 'waa_pointer_fleets', 'waa_pointer_indices'.")
			}
			if(is.null(WAA)) {
			  stop("'waa' (data) or 'WAA' (parameters) not provided, so there is no way to calculate population mean weight-at-age.")
			}
			n_pointers = max(c(waa_opts$waa_pointer_totcatch, waa_opts$waa_pointer_ssb, waa_opts$waa_pointer_jan1, waa_opts$waa_pointer_fleets, waa_opts$waa_pointer_indices)) 
			dim_WAA = c(n_pointers, data$n_years_model, data$n_ages)
					data$waa = array(2, dim = dim_WAA) # 2 kg for all ages, this will be replaced in WHAM	
					data$waa_cv = array(0, dim = dim_WAA)	# not use waa information for LL calculation	
		} else {
			data$waa = waa_opts$waa
			dim_waa = dim(data$waa)
			WAA_pointers = c(1:data$n_fleets,data$n_fleets+data$n_indices + c(1,2,2)) #Jan1 = SSB
			data$waa_pointer_indices = rep(1,data$n_indices) + data$n_fleets
			data$waa_pointer_fleets = WAA_pointers[1:data$n_fleets]
			data$waa_pointer_totcatch = WAA_pointers[data$n_fleets + 1]
			data$waa_pointer_ssb = WAA_pointers[data$n_fleets + 2]
			data$waa_pointer_jan1 = WAA_pointers[data$n_fleets + 3]
			data$waa_cv = array(0, dim = dim(data$waa))
		}
		
		if(!is.null(waa_opts$waa_pointer_fleets)) data$waa_pointer_fleets = waa_opts$waa_pointer_fleets
		if(!is.null(waa_opts$waa_pointer_indices)) data$waa_pointer_indices = waa_opts$waa_pointer_indices		
		if(!is.null(waa_opts$waa_pointer_totcatch)) data$waa_pointer_totcatch = waa_opts$waa_pointer_totcatch
		if(!is.null(waa_opts$waa_pointer_ssb)) data$waa_pointer_ssb = waa_opts$waa_pointer_ssb
		if(!is.null(waa_opts$waa_pointer_jan1)) data$waa_pointer_jan1 = waa_opts$waa_pointer_jan1
		if(!is.null(waa_opts$waa_cv)) data$waa_cv = waa_opts$waa_cv
		if(!is.null(waa_opts$use_catch_waa)) data$use_catch_waa = waa_opts$use_catch_waa
		if(!is.null(waa_opts$use_index_waa)) data$use_index_waa = waa_opts$use_index_waa

	}

  input$data = data

  # WAA parameters ----------------------------------------------------------------------------------

  data = input$data
  par = input$par
  map = input$map
    
  n_RE_cols = 4 # For RE parameters
  # sigma, rho_a, rho_y, rho_c
  
  # WAA default options:
  n_par_def = 2 
  # 2 parameters (LW): Omega1, Omega2
  # n_ages parameters (nonparametric): 1, ..., n_ages
  data$WAA_model = 1 # 1: LW, 2: Nonparametric
  data$WAA_est <- rep(0, times = n_par_def) # default = don't estimate WAA parameters
  data$n_WAA_par = n_par_def
  data$n_WAA_dim = 2  
  # 2 RE dim (LW): Omega1, Omega2
  # 1 RE dim (nonparametric)  
  data$WAA_re_model = rep(1, times = data$n_WAA_dim) # default = no RE / 'none'
  WAA_re_ini = array(0, dim = c(data$n_years_model, data$n_ages, data$n_WAA_dim))
  WAA_ini = c(log(5e-06), log(3))
  n_par = c(2, data$n_ages) # number of parameters for each method (LW: 2, nonparametric: n_ages)
  n_RE_dim = c(2, 1) # number of RE model for each method
  
  if(!is.null(WAA)) {

    # Use parametric option:
    data$isW_parametric = 1L
	
	# Some mandatory inputs should be provided:
	if(is.null(WAA$model)) stop("WAA$model must be provided.")
	if(is.null(WAA$init_vals)) stop("WAA$init_vals must be provided.")

	# WAA model:
    if(!(WAA$model %in% c("Allometric", "Nonparametric"))) stop("WAA$model must be 'Allometric' or 'Nonparametric'")
    data$WAA_model <- match(WAA$model, c("Allometric", "Nonparametric")) 
    data$n_WAA_par = n_par[data$WAA_model] # number of parameters to estimate
    data$WAA_est = rep(0, times = n_par[data$WAA_model]) # estimate?
    data$n_WAA_dim = n_RE_dim[data$WAA_model]
	data$WAA_re_model = rep(1, times = data$n_WAA_dim) # default = no RE / 'none'
    WAA_re_ini = array(0, dim = c(data$n_years_model, data$n_ages, n_RE_dim[data$WAA_model]))

	# WAA re model:
    if(!is.null(WAA$re)){
      if(length(WAA$re) != data$n_WAA_dim) stop("Length of 'WAA$re' must be equal to the number of RE dimensions of the chosen WAA model.")
      if(data$WAA_model %in% c(1)) {
		  for(k in 1:data$n_WAA_par) {
			  if(!(WAA$re[k] %in% c("none","iid","ar1"))) stop(paste0("WAA$re[", k, "] must be one of the following: 'none','iid','ar1'"))
			  data$WAA_re_model[k] <- match(WAA$re[k], c("none","iid","ar1")) # Respect this order to create array later
		  }
	  }
	  if(data$WAA_model %in% c(2)) {
		if(!(WAA$re %in% c("none","iid","2dar1","3dgmrf"))) stop("WAA$re must be one of the following: 'none','iid','2dar1','3dgmrf'")
		data$WAA_re_model <- match(WAA$re, c("none","iid","2dar1","3dgmrf")) # Respect this order to create array later			
	  }
    }
	
	# Init vals: mandatory
    if(length(WAA$init_vals) != data$n_WAA_par) stop("Length of 'WAA$init_vals' must be equal to the number of parameters of the chosen WAA model.")
    WAA_ini <- log(WAA$init_vals)

    if(!is.null(WAA$est_pars)){
        if(length(WAA$est_pars) > data$n_WAA_par) stop("Length of 'WAA$est_pars' should be equal or less than the number of parameters of the chosen WAA model.")
        data$WAA_est[WAA$est_pars] = 1
    }

  }

  # This section only used for 3D smoother:
  data$ay3D_IndexW = as.matrix(expand.grid("age" = seq_len(data$n_ages), "year" = seq_len(data$n_years_model) ))
  if(is.null(data$Var3D_ParamW)) data$Var3D_ParamW = 0 # Var_Param == 0 Conditional, == 1 Marginal

  # Update pars --------------------------
  
  par$WAA_a = WAA_ini
  par$WAA_re = WAA_re_ini
  par$WAA_repars = matrix(0, ncol = n_RE_cols, nrow = data$n_WAA_dim)
  par$WAA_repars[,1] = log(0.1) # start sigma at 0.1, rho at 0

  # --------------------------------
  # Map WAA:
  # Main pars:
  tmp1 <- par$WAA_a
  tmp1[data$WAA_est==0] = NA
  ind.notNA <- which(!is.na(tmp1))
  tmp1[ind.notNA] <- 1:length(ind.notNA)
  map$WAA_a <- factor(tmp1) 

  # RE info:
  map$WAA_re = NULL
  map$WAA_repars = NULL
  tmp.g.repars <- matrix(NA, ncol = n_RE_cols, nrow = data$n_WAA_dim)

  # Parametric equations:
  if(data$WAA_model %in% c(1)) {
    max_val_par = 0
	active_sum = -1
	this_max = 0
	for(i in 1:data$n_WAA_par) {
		tmp1 <- matrix(NA, nrow = data$n_years_model, ncol = data$n_ages)
		if(data$WAA_re_model[i] %in% c(2,3)){ # iid ar1- only y
		  tmp1[] = rep(1:nrow(tmp1), times = ncol(tmp1))  # all y estimated
		  active_sum = active_sum + 1
		  max_val_par = max_val_par + this_max * min(1, active_sum)
		  this_max = max(tmp1, na.rm = TRUE)
		}
		map$WAA_re <- c(map$WAA_re, as.vector(tmp1) + max_val_par)
		# WAA_repars: 
		if(data$WAA_re_model[i] == 1) tmp.g.repars[i,] <- rep(NA,n_RE_cols) # no RE pars to estimate
		if(data$WAA_re_model[i] == 2) tmp.g.repars[i,] <- c(1,NA,NA,NA) # estimate sigma 
		if(data$WAA_re_model[i] == 3) tmp.g.repars[i,] <- c(1,NA,2,NA) # estimate sigma ar1 y
	}
	map$WAA_re = as.factor(map$WAA_re)
	ind.notNA <- which(!is.na(tmp.g.repars))
	tmp.g.repars[ind.notNA] <- 1:length(ind.notNA)
	map$WAA_repars = factor(tmp.g.repars)
  }
  
  # Nonparametric:
  if(data$WAA_model %in% c(2)) { 
	  tmp <- par$WAA_re
	  if(data$WAA_re_model == 1) tmp[] = NA # no RE (either estimate RE for all ages or none at all)
	  if(data$WAA_re_model %in% 2:4){ # iid 2dar1 and 3dar1
		tmp[] = 1:(dim(tmp)[1]*dim(tmp)[2]) # all y,a estimated
	  }
	  map$WAA_re <- factor(tmp)
	  # WAA_repars: 
	  if(data$WAA_re_model == 1) tmp.g.repars <- rep(NA,n_RE_cols) # no RE pars to estimate
	  if(data$WAA_re_model == 2) tmp.g.repars <- c(1,NA,NA,NA) # estimate sigma
	  if(data$WAA_re_model == 3) tmp.g.repars <- c(1:3,NA) # 2dar1: estimate all
	  if(data$WAA_re_model == 4) tmp.g.repars <- c(1:4) # 3dar1: estimate all
	  map$WAA_repars = factor(tmp.g.repars)
  }

  # End section

  input$data = data
  input$par = par
  input$map = map

  return(input)
}
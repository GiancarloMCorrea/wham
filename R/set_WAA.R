set_WAA = function(input, waa_opts = NULL, WAA, LW)
{
  # Empirical Weight-at-age data ----------------------------------------------------------------------------------

  asap3 = if(is.null(input$asap3)) asap3 = NULL
  else asap3 = input$asap3
	data = input$data
  data$use_catch_waa = matrix(0, nrow = data$n_years_model, ncol = data$n_fleets)
  data$use_index_waa = matrix(0, nrow = data$n_years_model, ncol = data$n_indices)
  data$isW_ewaa = 1L # default value
	if(!is.null(asap3))
	{
	  i <- c(seq(1,(asap3$n_fleets+1)*2-1,2),(asap3$n_fleets+1)*2 + 1:2)
	  WAA_pointers <- asap3$WAA_pointers[i] #wham has no discard data, so remove those WAA matrices
	  data$waa = array(NA, dim = c(length(asap3$WAA_mats), data$n_years_model, data$n_ages))
	  for(i in 1:length(asap3$WAA_mats)) data$waa[i,,] = asap3$WAA_mats[[i]]
	  data$waa_pointer_indices = asap3$index_WAA_pointers
    data$waa_pointer_fleets = WAA_pointers[1:data$n_fleets]
    data$waa_pointer_totcatch = WAA_pointers[data$n_fleets + 1]
    data$waa_pointer_ssb = WAA_pointers[data$n_fleets + 2]
    data$waa_pointer_jan1 = WAA_pointers[data$n_fleets + 3]
    # if ASAP3 input provided, automatically use EWAA:
	  data$isW_ewaa = 1L
    data$isW_parametric = 0L
    data$isW_nonparametric = 0L
	  data$waa_cv = array(0, dim = dim(data$waa))
	} else {
    if(is.null(waa_opts[["waa"]])){
	    	#L = 100*(1-exp(-0.3*(1:data$n_ages - 0)))
		    #W = rep(exp(-11)*L^3, each = data$n_years_model)
        if(is.null(waa_opts$waa_pointer_totcatch) | is.null(waa_opts$waa_pointer_ssb) | is.null(waa_opts$waa_pointer_jan1) | is.null(waa_opts$waa_pointer_fleets) | is.null(waa_opts$waa_pointer_indices)) {
          stop("If 'waa' is not available, waa pointers must be provided: 
            'waa_pointer_totcatch', 'waa_pointer_ssb', 'waa_pointer_jan1', 'waa_pointer_fleets', 'waa_pointer_indices'.")
        }
        if(is.null(WAA) & is.null(LW) ) {
          stop("'waa' (data), 'LW', 'WAA' (parameters) not provided, so there is no way to calculate population mean weight-at-age. 
            Please, provide one of those.")
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
    
  n_re_par = 4 # number of parameters for RE

  # LAA default options:
  WAA_re_ini = matrix(0, ncol = data$n_ages, nrow = data$n_years_model)
  data$WAA_re_model = 1
  data$WAA_est = rep(0, times = data$n_ages)
  data$isW_nonparametric = 0L # default: it is not nonparametric
  WAA_ini = log( 2e-06*(100 + (3 - 100)*exp(-0.2*(1:data$n_ages - 1)))^3 )

  if(!is.null(LW)) data$isW_ewaa = 0L # turn off EWAA if LW provided

  if(!is.null(WAA)) {

    data$isW_nonparametric = 1L # it is nonparametric
    data$isW_ewaa = 0L # turn off EWAA if WAA provided

    WAA_ini = log(WAA$WAA_vals)
    if(!is.null(WAA$est_pars)) data$WAA_est[WAA$est_pars] = 1
    if(!is.null(WAA$re))  {
      if(!(WAA$re %in% c("none","iid","iid_a","ar1_a","2dar1","3dar1"))) stop("WAA$re must be one of the following: 'none','iid','iid_a','ar1_a','2dar1','3dar1'")
      data$WAA_re_model <- match(WAA$re, c("none","iid","iid_a","ar1_a","2dar1","3dar1")) # Respect this order to create array later
    }

  }

  # This section only used for 3D smoother:
  data$ay3D_IndexW = as.matrix(expand.grid("age" = seq_len(data$n_ages), "year" = seq_len(data$n_years_model) ))
  if(is.null(data$Var3D_ParamW)) data$Var3D_ParamW = 0 # Var_Param == 0 Conditional, == 1 Marginal

  # organize par --------------------------
  
  par$WAA_a = WAA_ini
  par$WAA_re = WAA_re_ini
  par$WAA_repars = rep(0, times = n_re_par)
  par$WAA_repars[1] = log(0.1) # start sigma at 0.1, rho at 0

  # --------------------------------
  # map WAA:
  map$WAA_a <- factor(rep(NA, times = data$n_ages))

  if(is.null(LW)) { # only do it if nonparametric
  
    tmp1 <- par$WAA_a
    tmp1[data$WAA_est==0] = NA
    ind.notNA <- which(!is.na(tmp1))
    tmp1[ind.notNA] <- 1:length(ind.notNA)
    map$WAA_a <- factor(tmp1)

  }

  # RE info:
  # LAA_re: "none","iid","iid_a","ar1_a","2dar1"
  tmp <- par$WAA_re
  if(data$WAA_re_model == 1) tmp[] = NA # no RE (either estimate RE for all ages or none at all)
  if(data$WAA_re_model %in% c(2,5,6)){ # 2d ar1
    tmp[] = 1:(dim(tmp)[1]*dim(tmp)[2]) # all y,a estimated
  }
  if(data$WAA_re_model %in% c(3,4)){ # ar1_a (devs by age, constant by year)
    for(i in 1:dim(tmp)[2]) tmp[,i] = i
  }
  map$WAA_re <- factor(tmp)

  # WAA_repars: 
  if(data$WAA_re_model == 1) tmp <- rep(NA,n_re_par) # no RE pars to estimate
  if(data$WAA_re_model == 2) tmp <- c(1,NA,NA, NA) # estimate sigma
  if(data$WAA_re_model == 3) tmp <- c(1,NA,NA,NA) # estimate sigma over ages
  if(data$WAA_re_model == 4) tmp <- c(1,2,NA,NA) # ar1_a: estimate sigma, rho_a
  if(data$WAA_re_model == 5) tmp <- c(1:3, NA) # 2dar1: estimate all
  if(data$WAA_re_model == 6) tmp <- 1:4 # 2dar1: estimate all
  map$WAA_repars = factor(tmp)

  # End section

  input$data = data
  input$par = par
  input$map = map

  return(input)
}
set_maturity = function(input, basic_info, maturity)
{
  data = input$data
  par = input$par
  map = input$map
  data$isMat_parametric = 0L # input data default

  # Maturity as data input --------------------------------------------------
  data$mature = t(matrix(1, data$n_ages, length(input$years)))
  if(!is.null(basic_info[["maturity"]])){
    if(!(length(basic_info[["maturity"]]) %in% c(1,data$n_ages*length(input$years)))) stop("basic_info$maturity has been specified, but its length is not 1 or length(ages)*length(years)")
    else data$mature[] = basic_info[["maturity"]]
  }
  
  # Maturity as parametric function -----------------------------------------
  
  n_re_par = 2 # number of parameters for RE

  # maturity default options:
  n_par_def = 2 # 2 parameters: Omega3 and Omega4: Mat = 1/(1+exp(-Omega3*(l-Omega4)))
  data$mat_re_model = rep(1, times = n_par_def) # default = no RE / 'none'
  data$n_mat_par = n_par_def # 
  data$mat_est <- rep(0, times = n_par_def) # default 
  mat_re_ini = array(0, dim = c(data$n_years_model, n_par_def))
  mat_ini = c(log(1), log(data$n_ages*0.5)) # age-logistic default (Omega3 and Omega4)
  data$mat_model = 1 # age-logistic default
  
  # prepare maturity options:
  if(!is.null(maturity)){
    
    # Use parametric option:
    data$isMat_parametric = 1L
    
    # Change initial values if len-specific:
    if(!is.null(maturity$model)) {
      if(maturity$model == 'len-logistic') { 
        mat_ini = c(log(2), log(max(data$lengths)*0.5)) # Length logistic default
        data$mat_model = 2
      }
    }

    if(!is.null(maturity$re)){
      
      if(length(maturity$re) != data$n_mat_par) stop("Number of 're' must be equal to the number of maturity parameters.")
      for(k in 1:data$n_mat_par) {
          if(!(maturity$re[k] %in% c("none","iid","ar1"))) stop(paste0("maturity$re[", k, "] must be one of the following: 'none','iid', or 'ar1'"))
          data$mat_re_model[k] <- match(maturity$re[k], c("none","iid","ar1")) # Respect this order to create array later
      }
      
    }

    if(!is.null(maturity$init_vals)){
      if(length(maturity$init_vals) != data$n_mat_par) stop(paste0("length(maturity$init_vals) must be ", n_par_def, "."))
      mat_ini <- log(maturity$init_vals)
    }
  
    if(!is.null(maturity$est_pars)){
        if(length(maturity$est_pars) > data$n_mat_par) stop(paste0("maturity$est_pars should contain values equal or less than ", n_par_def, "."))
        data$mat_est[maturity$est_pars] = 1
    }

  }
  data$n_mat_est <- sum(data$mat_est)

  # STOP if maturity as input and age-based maturity are provided
  if(!is.null(basic_info[["maturity"]]) & (data$isMat_parametric) & (data$mat_model == 1)) stop("Maturity as data input and age-logistic cannot be specified simultaneously.")

  # maturity pars --------------------------
  
  par$mat_a = as.matrix(mat_ini)
  par$mat_re = mat_re_ini
  par$mat_repars = matrix(0, ncol = n_re_par, nrow = data$n_mat_par)
  par$mat_repars[,1] = log(0.1) # start sigma at 0.1, rho at 0


  # --------------------------------
  # Prepare data for maturity:
  tmp1 <- par$mat_a

  tmp1[data$mat_est==0] = NA
  ind.notNA <- which(!is.na(tmp1))
  tmp1[ind.notNA] <- 1:length(ind.notNA)
  map$mat_a <- factor(tmp1)

  # RE info:
  map$mat_re = NULL
  map$mat_repars = NULL
  tmp.g.repars <- matrix(NA, ncol = n_re_par, nrow = data$n_mat_par)
  max_val_par = 0
  active_sum = -1
  this_max = 0
  for(i in 1:data$n_mat_par) {

    # mat_re: "none","iid","ar1"
    tmp1 <- rep(NA, times = data$n_years_model)
    if(data$mat_re_model[i] %in% c(2,3)){ # iid ar1- only y
      tmp1 = 1:length(tmp1) # all y estimated
      active_sum = active_sum + 1
      max_val_par = max_val_par + this_max * min(1, active_sum)
      this_max = max(tmp1, na.rm = TRUE)
    }

    map$mat_re <- c(map$mat_re, tmp1 + max_val_par)

    # K_repars: sigma_M, rho_M_y
    if(data$mat_re_model[i] == 1) tmp.g.repars[i,] <- rep(NA,n_re_par) # no RE pars to estimate
    if(data$mat_re_model[i] == 2) tmp.g.repars[i,] <- c(1,NA) # estimate sigma y
    if(data$mat_re_model[i] == 3) tmp.g.repars[i,] <- c(1,2) # estimate sigma and rho y

  }

  map$mat_re = as.factor(map$mat_re)
  ind.notNA <- which(!is.na(tmp.g.repars))
  tmp.g.repars[ind.notNA] <- 1:length(ind.notNA)
  map$mat_repars = factor(tmp.g.repars)

  # End section

  input$data = data
  input$par = par
  input$map = map
  
  return(input)

}
#define TMB_LIB_INIT R_init_wham
#include <TMB.hpp>
#include <iostream>
#include "helper_functions.hpp"
#include "age_comp_osa.hpp"
#include "age_comp_sim.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE

  DATA_INTEGER(n_years_catch); //same as n_years_model
  DATA_INTEGER(n_years_indices); //same as n_years_model
  DATA_INTEGER(n_years_model); 
  DATA_INTEGER(n_ages);
  DATA_INTEGER(n_lengths); 
  DATA_VECTOR(lengths); // uniform length bins
  DATA_INTEGER(n_fleets);
  DATA_INTEGER(n_indices);
  DATA_INTEGER(n_selblocks);
  DATA_IVECTOR(selblock_models); // for each block: 1 = age-specific, 2 = logistic, 3 = double-logistic, 4 = logistic (declining), 5 = age double normal, 6 = len logistic, 7 = len decreasing logistic, 8 = len double normal
  DATA_IVECTOR(selblock_models_re); // for each block: 1 = none, 2 = IID, 3 = ar1, 4 = ar1_y, 5 = 2dar1
  DATA_IVECTOR(n_selpars);
  DATA_IMATRIX(selpars_est); // n_blocks x (n_pars + n_ages), is the selpar estimated in this block?
  DATA_IVECTOR(n_selpars_est); // of the selpars, how many are actually estimated (not fixed at 0 or 1)
  DATA_IVECTOR(n_years_selblocks); // for each block, number of years the block covers
  DATA_IMATRIX(selblock_years); // n_years_model x n_selblocks, = 1 if block covers year, = 0 if not
  DATA_IMATRIX(selblock_pointer_fleets);
  DATA_IMATRIX(selblock_pointer_indices);
  DATA_IVECTOR(age_comp_model_fleets);
  DATA_IVECTOR(age_comp_model_indices);
  DATA_IVECTOR(len_comp_model_fleets); 
  DATA_IVECTOR(len_comp_model_indices); 
  DATA_VECTOR(fracyr_SSB);
  DATA_MATRIX(mature);
  DATA_IVECTOR(waa_pointer_fleets);
  DATA_INTEGER(waa_pointer_totcatch);
  DATA_IVECTOR(waa_pointer_indices);
  DATA_INTEGER(waa_pointer_ssb);
  DATA_INTEGER(waa_pointer_jan1);
  DATA_IMATRIX(use_catch_waa);
  DATA_IMATRIX(use_index_waa);
  DATA_ARRAY(waa);
  DATA_ARRAY(waa_cv);
  DATA_MATRIX(agg_catch);
  DATA_VECTOR(fracyr_catch);
  DATA_IMATRIX(use_agg_catch);
  DATA_MATRIX(agg_catch_sigma);
  DATA_ARRAY(catch_paa); //n_fleets x n_years x n_ages
  DATA_IMATRIX(use_catch_paa);
  DATA_MATRIX(catch_Neff);
  DATA_ARRAY(catch_aging_error); 
  DATA_VECTOR(use_catch_aging_error); 
  DATA_ARRAY(catch_pal); //n_fleets x n_years x n_lengths 
  DATA_IMATRIX(use_catch_pal); 
  DATA_MATRIX(catch_NeffL); 
  DATA_ARRAY(catch_caal); //n_fleets x n_years x n_lengths x n_ages
  DATA_IARRAY(use_catch_caal); 
  DATA_ARRAY(catch_caal_Neff);
  DATA_IVECTOR(units_indices);
  DATA_MATRIX(fracyr_indices);
  DATA_MATRIX(agg_indices);
  DATA_IMATRIX(use_indices);
  DATA_MATRIX(agg_index_sigma);
  DATA_IVECTOR(units_index_paa);
  DATA_ARRAY(index_paa); //n_indices x n_years x n_ages
  DATA_IMATRIX(use_index_paa);
  DATA_MATRIX(index_Neff);
  DATA_IVECTOR(units_index_pal); 
  DATA_ARRAY(index_pal); //n_indices x n_years x n_lengths 
  DATA_IMATRIX(use_index_pal); 
  DATA_MATRIX(index_NeffL);
  DATA_ARRAY(index_caal); //n_fleets x n_years x n_lengths x n_ages
  DATA_IARRAY(use_index_caal); 
  DATA_ARRAY(index_caal_Neff); 
  DATA_ARRAY(index_aging_error);   
  DATA_VECTOR(use_index_aging_error); 
  DATA_VECTOR(q_lower);
  DATA_VECTOR(q_upper);
  DATA_IVECTOR(use_q_prior);
  DATA_VECTOR(logit_q_prior_sigma);
  DATA_IVECTOR(use_q_re); //n_indices, 0= no re, >0 = use re 
  DATA_MATRIX(selpars_lower);
  DATA_MATRIX(selpars_upper);
  DATA_INTEGER(n_NAA_sigma); // 0 = SCAA, 1 = logR only, 2 = full state-space with shared sig_a for a > 1
  DATA_IVECTOR(NAA_sigma_pointers);
  DATA_INTEGER(recruit_model);
  DATA_INTEGER(n_M_a);
  DATA_INTEGER(M_model); // 1: "constant", 2: "age-specific", 3: "weight-at-age"
  DATA_INTEGER(N1_model); //0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations
  DATA_INTEGER(M_re_model); // 1 = none, 2 = IID, 3 = ar1_a, 4 = ar1_y, 5 = 2dar1
  // M_est and n_M_est are not use
  //DATA_IVECTOR(M_est); // Is mean M estimated for each age? If M-at-age, dim = length(n_M_a). If constant or weight-at-age M, dim = 1.
  //DATA_INTEGER(n_M_est); // How many mean M pars are estimated?
  DATA_INTEGER(use_b_prior);
  // WAA information
  DATA_INTEGER(n_WAA_dim); // Number of RE dimensions  
  DATA_INTEGER(WAA_model); // 1: "Allometric", 2: "Nonparametric"
  DATA_IVECTOR(WAA_re_model); // Depends on the WAA model
  DATA_IVECTOR(WAA_est); 
  DATA_INTEGER(isW_parametric); // WAA model is parametric? 1 = yes, 0 = no
  DATA_MATRIX(ay3D_IndexW);  // (n_years * n_ages) * 2 
  DATA_INTEGER(Var3D_ParamW); // Variance parameterization of Precision Matrix == 0 (Conditional), == 1(Marginal)
  // LAA information
  DATA_INTEGER(n_LAA_par); // Number of fixed effects
  DATA_INTEGER(n_LAA_dim); // Number of RE dimensions  
  DATA_INTEGER(LAA_model); // 1: "vB-classic", 2: "Richards", 3: "Nonparametric"
  DATA_IVECTOR(LAA_re_model); // Depends on the LAA model
  DATA_IVECTOR(LAA_est); 
  DATA_SCALAR(age_L1); // age for L1
  DATA_INTEGER(age_L1_ceil); // age (ceiling) for L1
  DATA_MATRIX(ay3D_IndexL);  // (n_years * n_ages) * 2 
  DATA_INTEGER(Var3D_ParamL); // Variance parameterization of Precision Matrix == 0 (Conditional), == 1(Marginal)
  // Maturity information
  DATA_INTEGER(isMat_parametric); // Maturity model is parametric? 1 = yes, 0 = no  
  DATA_INTEGER(mat_model); // age or length-specific
  DATA_INTEGER(n_mat_par);
  DATA_IVECTOR(mat_re_model); // 1 = none, 2 = IID_y, 3 = ar1_y
  // Continue..
  DATA_IVECTOR(which_F_age); //which age of F to use for full total F for msy/ypr calculations and projections (n_years_model + n_years_proj)
  DATA_INTEGER(use_steepness); // which parameterization to use for BH/Ricker S-R, if needed.
  DATA_INTEGER(bias_correct_pe); //bias correct lognormal process error?
  DATA_INTEGER(bias_correct_oe); //bias correct lognormal observation error?
  DATA_IVECTOR(Fbar_ages);
  DATA_IVECTOR(simulate_state); //vector (0/1) if 1 then state parameters (NAA, MAA, sel, Ecov, q, LAA, WAA, Maturity) in that order) will be simulated.
  DATA_IVECTOR(simulate_data); //vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
  DATA_IVECTOR(simulate_period); //vector (0/1) if 1 then period (model years, projection years) will be simulated.
  DATA_SCALAR(percentSPR); // percentage to use for SPR-based reference points. Default = 40.
  DATA_SCALAR(percentFXSPR); // percent of F_XSPR to use for calculating catch in projections. For example, GOM cod uses F = 75% F_40%SPR, so percentFXSPR = 75 and percentSPR = 40. Default = 100.
  DATA_SCALAR(percentFMSY); // percent of FMSY to use for calculating catch in projections.
  DATA_INTEGER(XSPR_R_opt); //1(3): use annual R estimates(predictions) for annual SSB_XSPR, 2(4): use average R estimates(predictions). See next line for years to average over.
  DATA_IVECTOR(XSPR_R_avg_yrs); // model year indices (TMB, starts @ 0) to use for averaging recruitment when defining SSB_XSPR (if XSPR_R_opt = 2,4)
  DATA_VECTOR(FXSPR_init); // annual initial values to use for newton steps to find FXSPR (n_years_model+n_proj_years)
  DATA_VECTOR(FMSY_init); // annual initial values to use for newton steps to find FMSY (n_years_model+n_proj_years)

  // data for one-step-ahead (OSA) residuals
  DATA_INTEGER(do_osa); //whether to do osa residuals. For efficiency reasons with age comp likelihoods.
  DATA_VECTOR(obsvec); // vector of all observations for OSA residuals
  DATA_IVECTOR(agesvec);
  DATA_IVECTOR(lensvec);																													  
  DATA_VECTOR_INDICATOR(keep, obsvec); // for OSA residuals
  DATA_IMATRIX(keep_C); // indices for catch obs, can loop years/fleets with keep(keep_C(y,f))
  DATA_IMATRIX(keep_I);
  DATA_IMATRIX(keep_E); // Ecov
  DATA_IARRAY(keep_Cpaa);
  DATA_IARRAY(keep_Ipaa);
  DATA_IARRAY(keep_Cwaa);
  DATA_IARRAY(keep_Iwaa);
  DATA_IARRAY(keep_Cpal); 
  DATA_IARRAY(keep_Ipal); 
  DATA_IARRAY(keep_Ccaal); 
  DATA_IARRAY(keep_Icaal); 
  DATA_IVECTOR(do_post_samp); //length = 9, whether to ADREPORT posterior residuals for NAA, M, selectivity, Ecov, q, LAA, WAA, Maturity

  // data for environmental covariate(s), Ecov
  DATA_INTEGER(n_Ecov); // also = 1 if no Ecov
  DATA_INTEGER(n_years_Ecov); // num years in Ecov  process model
  DATA_IMATRIX(Ecov_use_obs); // all 0 if no Ecov
  DATA_MATRIX(Ecov_obs);
  //Below is not used anymore
  //DATA_IVECTOR(Ecov_lag);
  DATA_IVECTOR(Ecov_how); // specific to recruitment effects. 0 = no effect, 1 = controlling, 2 = limiting, 3 = lethal, 4 = masking, 5 = directive
  //Below is not used anymore
  //DATA_IMATRIX(Ecov_poly); // n_Ecov x 2+n_indices. polynomial order for ecov effects (1 = linear, 2 = quadratic, 3 = cubic, ...)
  DATA_IMATRIX(Ecov_where); // n_Ecov x 2+n_indices. 0/1 values with columns corresponding to recruit, mortality, indices, LAA, WAA, Maturity in that order
  DATA_IVECTOR(Ecov_model); // 0 = no Ecov, 1 = RW, 2 = AR1
  //DATA_INTEGER(year1_Ecov); // first year Ecov
  //DATA_INTEGER(year1_model); // first year model
  DATA_IMATRIX(ind_Ecov_out_start); // n_Ecov x (2 + n_indices) index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects specific the multiple types of effects each Ecov can have)
  DATA_IMATRIX(ind_Ecov_out_end); // n_Ecov x (2 + n_indices) index of Ecov_x to use for Ecov_out (operates on pop model, lagged effects specific the multiple types of effects each Ecov can have)
  DATA_IVECTOR(Ecov_obs_sigma_opt); // n_Ecov, 1 = given, 2 = estimate 1 value, shared among obs, 3 = estimate for each obs, 4 = estimate for each obs as random effects
  DATA_IMATRIX(Ecov_use_re); // 0/1: use Ecov_re? If yes, add to nll. (n_years_Ecov + n_years_proj_Ecov) x n_Ecov

  // data for projections
  DATA_INTEGER(do_proj); // 1 = yes, 0 = no
  DATA_INTEGER(n_years_proj); // number of years to project
  DATA_INTEGER(n_years_proj_Ecov); // number of years to project Ecov
  DATA_IVECTOR(avg_years_ind); // model year indices (TMB, starts @ 0) to use for averaging MAA, waa, maturity, and F (if use.avgF = TRUE)
  DATA_IVECTOR(proj_F_opt); // for each projection year, 1 = last year F (default), 2 = average F, 3 = F at X% SPR, 4 = user-specified F, 5 = calculate F from user-specified catch
  DATA_VECTOR(proj_Fcatch); // user-specified F or catch in projection years, only used if proj_F_opt = 4 or 5
  DATA_INTEGER(proj_M_opt); // 1 = continue M_re (check for time-varying M_re on R side), 2 = average M (over avg_years_ind)
  DATA_IVECTOR(proj_LAA_opt); // 1 = continue LAA_re (check for time-varying LAA_re on R side), 2 = average LAA (over avg_years_ind)
  DATA_IVECTOR(proj_WAA_opt); // 1 = continue WAA_re (check for time-varying WAA_re on R side), 2 = average WAA (over avg_years_ind)
  DATA_IVECTOR(proj_mat_opt); // 1 = continue mat_re (check for time-varying mat_re on R side), 2 = average mat (over avg_years_ind)
  DATA_SCALAR(logR_mean); // empirical mean recruitment in model years, used for SCAA recruit projections
  DATA_SCALAR(logR_sd); // empirical sd recruitment in model years, used for SCAA recruit projections
  DATA_VECTOR(F_proj_init); // annual initial values  to use for newton steps to find F for use in projections  (n_years_proj)
  
  //static brp info
  DATA_INTEGER(which_F_age_static); //which age of F to use for full total F for static brps (max of average FAA_tot over avg_years_ind)
  DATA_SCALAR(static_FXSPR_init); // initial value to use for newton steps to find FXSPR_static
  //DATA_IVECTOR(static_brp_years_sel); // model year indices (TMB, starts @ 0) to use for averaging selectivity for static biological reference points
  //DATA_IVECTOR(static_brp_years_mat); // model year indices (TMB, starts @ 0) to use for averaging maturity for static biological reference points
  //DATA_IVECTOR(static_brp_years_waa_ssb); // model year indices (TMB, starts @ 0) to use for averaging SSB weight at age for static biological reference points
  //DATA_IVECTOR(static_brp_years_waa_catch); // model year indices (TMB, starts @ 0) to use for averaging catch weight at age for static biological reference points
  //DATA_IVECTOR(static_brp_years_M); // model year indices (TMB, starts @ 0) to use for averaging nat mort. for static biological reference points
  //DATA_IVECTOR(static_brp_years_R); // ******just use XSPR_R_avg_yrs********model year indices (TMB, starts @ 0) to use for averaging recruitment for static biological reference points

  // parameters - general
  PARAMETER_VECTOR(mean_rec_pars);
  PARAMETER_VECTOR(logit_q); //n_indices (mean/constant q pars)
  PARAMETER_VECTOR(q_prior_re); //n_indices (if a prior is used for q, this is used instead of logit_q)
  PARAMETER_ARRAY(q_re); //n_years x n_indices (time series of)
  PARAMETER_MATRIX(q_repars) //n_indices x 2 (sigma, rho)
  PARAMETER_VECTOR(log_F1);
  PARAMETER_MATRIX(F_devs);
  PARAMETER_VECTOR(log_N1_pars); //length = n_ages or 2
  PARAMETER_VECTOR(log_NAA_sigma);
  PARAMETER_VECTOR(trans_NAA_rho); // rho_a, rho_y (length = 2)
  PARAMETER_ARRAY(log_NAA); // numbers-at-age random effects, (dim = n.yrs-1, n.ages)
  PARAMETER_VECTOR(logR_proj); // recruitment (random effects) in proj years, only if SCAA
  // PARAMETER_ARRAY(NAA_re); // random effects / deviations from pred NAA for year y, age a (dim = n.yrs-1, n.ages)
  PARAMETER_MATRIX(logit_selpars); // mean selectivity, dim = n_selblocks x n_ages + 6 (n_ages for by age, 2 for logistic, 4 for double-logistic)
  PARAMETER_VECTOR(selpars_re);    // deviations in selectivity parameters (random effects), length = sum(n_selpars)*n_years per block
  PARAMETER_MATRIX(sel_repars);    // fixed effect parameters controlling selpars_re, dim = n_blocks, 3 (sigma, rho, rho_y)
  //PARAMETER_VECTOR(catch_paa_pars);
  PARAMETER_MATRIX(catch_paa_pars); //n_fleets x 3
  PARAMETER_MATRIX(catch_pal_pars); //n_fleets x 3 
  //PARAMETER_VECTOR(index_paa_pars);
  PARAMETER_MATRIX(index_paa_pars); //n_indices x 3
  PARAMETER_MATRIX(index_pal_pars); //n_indices x 3 
  PARAMETER_VECTOR(M_a); // mean M-at-age, fixed effects, length = n_ages if M_model = 2 (age-specific), length = 1 if M_model = 1 (constant) or 3 (weight-at-age M)
  PARAMETER_ARRAY(M_re); // random effects for year- and age-varying M deviations from mean M_a), dim = n_years x n_M_a
  PARAMETER_VECTOR(M_repars); // parameters controlling M_re, length = 3 (sigma_M, rho_M_a, rho_M_y)
  // LAA parameters:
  PARAMETER_VECTOR(LAA_a);  
  PARAMETER_ARRAY(LAA_re);    
  PARAMETER_MATRIX(LAA_repars);    
  PARAMETER_VECTOR(SD_par);
  // WAA parameters:
  PARAMETER_VECTOR(WAA_a);  
  PARAMETER_ARRAY(WAA_re);  
  PARAMETER_MATRIX(WAA_repars);  
  // Maturity parameters:  
  PARAMETER_VECTOR(mat_a);  
  PARAMETER_ARRAY(mat_re); 
  PARAMETER_MATRIX(mat_repars);  
  // Continue..   
  PARAMETER(log_b);
  PARAMETER_VECTOR(log_catch_sig_scale); //n_fleets
  PARAMETER_VECTOR(log_index_sig_scale); //n_indices

  // parameters - environmental covariate ("Ecov")
  PARAMETER_MATRIX(Ecov_re); // nrows = n_years_Ecov, ncol = N_Ecov
  PARAMETER_ARRAY(Ecov_beta); // dim = (2 + n_indices) x n_poly x n_ecov x n_ages, beta_R in eqns 4-5, Miller et al. (2016)
  PARAMETER_MATRIX(Ecov_process_pars); // nrows = RW: 2 par (Ecov1, sig), AR1: 3 par (mu, sig, phi); ncol = N_ecov
  PARAMETER_MATRIX(Ecov_obs_logsigma); // N_Ecov_years x n_Ecov. options: just given (data), or fixed effect(s)
  PARAMETER_MATRIX(Ecov_obs_logsigma_re); // N_Ecov_years x n_Ecov. columns of random effects used if Ecov_obs_sigma_opt = 4 
  PARAMETER_MATRIX(Ecov_obs_sigma_par); // ncol = N_Ecov, nrows = 2 (mean, sigma of random effects)

  Type nll= 0.0; //negative log-likelihood
  vector<int> any_index_age_comp(n_indices);
  vector<int> any_fleet_age_comp(n_fleets);
  vector<int> any_index_len_comp(n_indices);  
  vector<int> any_fleet_len_comp(n_fleets);  
  vector<int> any_index_caal(n_indices);  
  vector<int> any_fleet_caal(n_fleets);  
  vector<Type> SSB(n_years_model + n_years_proj);
  matrix<Type> F(n_years_model,n_fleets);
  matrix<Type> log_F(n_years_model,n_fleets);
  
  array<Type> pred_CAA(n_years_model+n_years_proj,n_fleets,n_ages);
  array<Type> pred_catch_paa(n_years_model+n_years_proj,n_fleets,n_ages);
  array<Type> pred_CAL(n_years_model+n_years_proj,n_fleets,n_lengths); 
  array<Type> pred_CAAL(n_years_model+n_years_proj,n_fleets,n_lengths,n_ages);
  array<Type> pred_catch_pal(n_years_model+n_years_proj,n_fleets,n_lengths); 
  array<Type> pred_catch_caal(n_years_model+n_years_proj,n_fleets,n_lengths,n_ages); 
  matrix<Type> pred_catch(n_years_model+n_years_proj,n_fleets);
  matrix<Type> pred_log_catch(n_years_model+n_years_proj,n_fleets);
  
  array<Type> pred_IAA(n_years_model+n_years_proj,n_indices,n_ages);
  array<Type> pred_IAAL(n_years_model+n_years_proj,n_indices,n_lengths,n_ages); 
  array<Type> pred_index_paa(n_years_model+n_years_proj,n_indices,n_ages);
  array<Type> pred_IAL(n_years_model+n_years_proj,n_indices,n_lengths);
  array<Type> pred_index_pal(n_years_model+n_years_proj,n_indices,n_lengths); 
  array<Type> pred_index_caal(n_years_model+n_years_proj,n_indices,n_lengths,n_ages); 
  matrix<Type> pred_indices(n_years_model+n_years_proj,n_indices); // not bias corrected
  matrix<Type> pred_log_indices(n_years_model+n_years_proj,n_indices); // bias corrected  

  array<Type> pred_waa(waa.dim(0), n_years_model + n_years_proj, n_ages); 
  array<Type> jan1_phi_mat(n_lengths,n_ages,n_years_model + n_years_proj);
  matrix<Type> jan1LAA(n_years_model + n_years_proj,n_ages); // mean length-at-age on Jan 1st
  matrix<Type> SDAA(n_years_model + n_years_proj,n_ages);
  vector<Type> temp_selAA(n_ages);
  array<Type> FAA(n_years_model+n_years_proj,n_fleets,n_ages);
  array<Type> log_FAA(n_years_model+n_years_proj,n_fleets,n_ages);
  matrix<Type> FAA_tot(n_years_model + n_years_proj,n_ages);
  matrix<Type> ZAA(n_years_model + n_years_proj,n_ages);
  // array<Type> QAA(n_years_model+n_years_proj,n_indices,n_ages);
  vector<array<Type> > phi_matrix(waa.dim(0)); // save phi matrix at different fracyr. Each array = y,l,a
  vector<matrix<Type> > selAL(n_selblocks); // Could be either selex-at-age or selex-at-len
  vector<matrix<Type> > selAA(n_selblocks); // selAA(b)(y,a) gives selectivity by block, year, age; selAA(b) is matrix with dim = n_years x n_ages;
  vector<matrix<Type> > selLL(n_selblocks); // Only selex-at-len. Put 1 for all length bins when main selex is selex-at-age
  matrix<Type> q(n_years_model+n_years_proj,n_indices);
  vector<Type> t_paa(n_ages); 
  vector<Type> t_pred_paa(n_ages); 
  vector<Type> t_pred_pal(n_lengths); 
  int n_toavg = avg_years_ind.size();

  //Type SR_a, SR_b, SR_R0, SR_h;
  for(int i = 0; i < n_indices; i++)
  {
    any_index_age_comp(i) = 0;
    for(int y = 0; y < n_years_indices; y++) if(use_index_paa(y,i) == 1) any_index_age_comp(i) = 1;
  }
  for(int i = 0; i < n_fleets; i++)
  {
    any_fleet_age_comp(i) = 0;
    for(int y = 0; y < n_years_catch; y++) if(use_catch_paa(y,i) == 1) any_fleet_age_comp(i) = 1;
  }
  for(int i = 0; i < n_indices; i++)
  {
    any_index_len_comp(i) = 0;
    for(int y = 0; y < n_years_indices; y++) if(use_index_pal(y,i) == 1) any_index_len_comp(i) = 1;
  }
  for(int i = 0; i < n_fleets; i++)
  {
    any_fleet_len_comp(i) = 0;
    for(int y = 0; y < n_years_catch; y++) if(use_catch_pal(y,i) == 1) any_fleet_len_comp(i) = 1;
  }
  for(int i = 0; i < n_indices; i++)
  {
    any_index_caal(i) = 0;
    for(int y = 0; y < n_years_indices; y++) {
		for(int l = 0; l < n_lengths; l++) if(use_index_caal(y,i,l) == 1) any_index_caal(i) = 1;
	}
  }
  for(int f = 0; f < n_fleets; f++)
  {
    any_fleet_caal(f) = 0;
    for(int y = 0; y < n_years_catch; y++) {
		for(int l = 0; l < n_lengths; l++) if(use_catch_caal(y,f,l) == 1) any_fleet_caal(f) = 1;
	}
  }

  // Selectivity --------------------------------------------------------------
  vector<array<Type> > selpars_re_mats(n_selblocks); // gets selectivity deviations (RE vector, selpars_re) as vector of matrices (nyears x npars), one for each block
  vector<matrix<Type> > selpars(n_selblocks); // selectivity parameter matrices for each block, nyears x npars
  Type nll_sel = Type(0);
  int istart = 0;
  int ct = 0;
  for(int b = 0; b < n_selblocks; b++){
    array<Type> tmp2(n_years_model,n_selpars(b));
    tmp2.setZero();
    selpars_re_mats(b) = tmp2;

    int jstart = 0; // offset for indexing selectivity pars, depends on selectivity model for block b: n_ages (age-specific) + 2 (logistic) + 4 (double-logistic)
    if((selblock_models(b) == 2) | (selblock_models(b) == 4)) jstart = n_ages;
    if(selblock_models(b) == 3) jstart = n_ages + 2; // 
    if(selblock_models(b) == 5) jstart = n_ages + 6;
    if((selblock_models(b) == 6) | (selblock_models(b) == 7)) jstart = n_ages + 12;
    if(selblock_models(b) == 8) jstart = n_ages + 14;

    if(selblock_models_re(b) > 1){
      // fill in sel devs from RE vector, selpars_re (fixed at 0 if RE off)
      array<Type> tmp(n_years_selblocks(b), n_selpars_est(b));
      for(int j=0; j<n_selpars_est(b); j++){
        tmp.col(j) = selpars_re.segment(istart,n_years_selblocks(b));
        istart += n_years_selblocks(b);
      }

      //question: is it faster here to just work on the selectivity parameters as re rather than the deviations?
      // likelihood of RE sel devs (if turned on)
      Type sigma; // sd selectivity deviations (fixed effect)
      Type rho; // among-par correlation selectivity deviations (fixed effect)
      Type rho_y; // among-year correlation selectivity deviations (fixed effect)
      Type Sigma_sig_sel;
      sigma = exp(sel_repars(b,0));
      rho = rho_trans(sel_repars(b,1)); // rho_trans ensures correlation parameter is between -1 and 1, see helper_functions.hpp
      rho_y = rho_trans(sel_repars(b,2)); // rho_trans ensures correlation parameter is between -1 and 1, see helper_functions.hpp
      
      if((selblock_models_re(b) == 2) | (selblock_models_re(b) == 5)){
        // 2D AR1 process on selectivity parameter deviations
        Sigma_sig_sel = pow(pow(sigma,2) / ((1-pow(rho_y,2))*(1-pow(rho,2))),0.5);
        nll_sel += SCALE(SEPARABLE(AR1(rho),AR1(rho_y)), Sigma_sig_sel)(tmp);
        SIMULATE if(simulate_state(2) == 1) if(sum(simulate_period) > 0) SEPARABLE(AR1(rho),AR1(rho_y)).simulate(tmp);
      } else {
        // 1D AR1 process on selectivity parameter deviations
        if(selblock_models_re(b) == 3){ // ar1 across parameters in selblock, useful for age-specific pars.
          vector<Type> tmp0 = tmp.matrix().row(0); //random effects are constant across years 
          Sigma_sig_sel = pow(pow(sigma,2) / (1-pow(rho,2)),0.5);
          nll_sel += SCALE(AR1(rho), Sigma_sig_sel)(tmp0);
          SIMULATE if(simulate_state(2) == 1) if(sum(simulate_period) > 0) 
          {
            AR1(rho).simulate(tmp0);
            for(int y = 0; y < tmp.rows(); y++) for(int i = 0; i < tmp0.size(); i++) tmp(y,i) = tmp0(i);
          }
        } else { // selblock_models_re(b) = 4, ar1_y, not sure if this one really makes sense.
          vector<Type> tmp0 = tmp.matrix().col(0); //random effects are constant within years 
          Sigma_sig_sel = pow(pow(sigma,2) / (1-pow(rho_y,2)),0.5);
          //Sigma_sig_sel = sigma;
          nll_sel += SCALE(AR1(rho_y), Sigma_sig_sel)(tmp0);
          SIMULATE if(simulate_state(2) == 1) if(sum(simulate_period) > 0)  
          {
            AR1(rho_y).simulate(tmp0);
            for(int a = 0; a < tmp.cols(); a++) tmp.col(a) = tmp0;
          }
        }
      }
      SIMULATE if(simulate_state(2) == 1) if(sum(simulate_period) > 0) {
        tmp = tmp * Sigma_sig_sel;
        istart -= n_selpars_est(b) * n_years_selblocks(b); //bring it back to the beginning for this selblock
        for(int j=0; j<n_selpars_est(b); j++){
          for(int y = 0; y < n_years_selblocks(b); y++){
            selpars_re(istart) = tmp(y,j);
            istart++;
          }
        }
      }

      // construct deviations array with full dimensions (n_years_model instead of n_years_selblocks, n_selpars instead of n_selpars_est)
      for(int j=0; j<n_selpars(b); j++){
        for(int y=0; y<n_years_model; y++){
          if((selblock_years(y,b) == 1) & (selpars_est(b,j+jstart) > 0)){
            selpars_re_mats(b)(y,j) = selpars_re(ct);
            ct++;
          }
        }
      }
    }

    // get selpars = mean + deviations
    matrix<Type> tmp1(n_years_model, n_selpars(b));
    for(int j=jstart; j<(jstart+n_selpars(b)); j++){ // transform from logit-scale
      for(int i=0; i<n_years_model; i++){
        tmp1(i,j-jstart) = selpars_lower(b,j) + (selpars_upper(b,j) - selpars_lower(b,j)) / (1.0 + exp(-(logit_selpars(b,j) + selpars_re_mats(b).matrix()(i,j-jstart))));
      }
    }
    selpars(b) = tmp1;
  }
  REPORT(selpars);
  REPORT(sel_repars);
  REPORT(selpars_re); //even if not simulated
  if(do_post_samp(2) == 1) ADREPORT(selpars_re);
  REPORT(logit_selpars);
  REPORT(nll_sel);
  selAL = get_selectivity(n_years_model, n_ages, n_lengths, lengths, n_selblocks, selpars, selblock_models); // Get selectivity by block, age, year. This contains either selex-at-age or selex-at-len
  nll += nll_sel;

  // Environmental covariate process model --------------------------------------
  matrix<Type> Ecov_x(n_years_Ecov + n_years_proj_Ecov, n_Ecov); // 'true' estimated Ecov (x_t in Miller et al. 2016 CJFAS)
  matrix<Type> nll_Ecov(n_years_Ecov + n_years_proj_Ecov, n_Ecov); // nll contribution each Ecov_re
  nll_Ecov.setZero();

  // Ecov_model == 0) no Ecov
  for(int i = 0; i < n_Ecov; i++){ // loop over Ecovs
    // Ecov model option 1: RW
    if(Ecov_model(i) == 1){
      Type Ecov_sig; // sd (sig_x in Eq1, pg 1262, Miller et al. 2016)
      Ecov_sig = exp(Ecov_process_pars(1,i));
      Type Ecov1; // Ecov_x in year 1 (fixed effect)
      Ecov1 = Ecov_process_pars(0,i);

      Ecov_x(0,i) = Ecov1;
      nll_Ecov(1,i) -= dnorm(Ecov_re(1,i), Ecov1, Ecov_sig, 1); // Ecov_re(0,i) set to NA
      SIMULATE if((simulate_state(3) == 1) & (Ecov_use_re(1,i) == 1)) {
        if(simulate_period(0) == 1) {
          Ecov_re(1,i) = rnorm(Ecov1, Ecov_sig);
        }
      }
      Ecov_x(1,i) = Ecov_re(1,i);
      for(int y = 2; y < n_years_Ecov + n_years_proj_Ecov; y++){
        nll_Ecov(y,i) -= dnorm(Ecov_re(y,i), Ecov_re(y-1,i), Ecov_sig, 1);
        SIMULATE if((simulate_state(3) == 1) & (Ecov_use_re(y,i) == 1)) {
          if(((simulate_period(0) == 1) & (y < n_years_Ecov)) | ((simulate_period(1) == 1) & (y > n_years_Ecov-1))) {
            Ecov_re(y,i) = rnorm(Ecov_re(y-1,i), Ecov_sig);
          }
        }
        Ecov_x(y,i) = Ecov_re(y,i);
      }
    }

    // Ecov model option 2: AR1
    if(Ecov_model(i) == 2){
      Type Ecov_mu; // mean
      Type Ecov_phi; // autocorrelation
      Type Ecov_sig; // conditional sd
      Ecov_mu = Ecov_process_pars(0,i);
      Ecov_phi = -Type(1) + Type(2)/(Type(1) + exp(-Ecov_process_pars(2,i)));
      Ecov_sig = exp(Ecov_process_pars(1,i));

      nll_Ecov(0,i) -= dnorm(Ecov_re(0,i), Type(0), Ecov_sig*exp(-Type(0.5) * log(Type(1) - pow(Ecov_phi,Type(2)))), 1);
      SIMULATE if((simulate_state(3) == 1) & (Ecov_use_re(0,i) == 1)) {
        if(simulate_period(0) == 1) {
          Ecov_re(0,i) = rnorm(Type(0), Ecov_sig*exp(-Type(0.5) * log(Type(1) - pow(Ecov_phi,Type(2)))));
        }
      }
      for(int y = 1; y < n_years_Ecov + n_years_proj_Ecov; y++)
      {
        nll_Ecov(y,i) -= dnorm(Ecov_re(y,i), Ecov_phi * Ecov_re(y-1,i), Ecov_sig, 1);
        SIMULATE if((simulate_state(3) == 1) & (Ecov_use_re(y,i) == 1)) {
          if(((simulate_period(0) == 1) & (y < n_years_Ecov)) | ((simulate_period(1) == 1) & (y > n_years_Ecov-1))) {
            Ecov_re(y,i) = rnorm(Ecov_phi * Ecov_re(y-1,i), Ecov_sig);
          }
        }
      }
      for(int y = 0; y < n_years_Ecov + n_years_proj_Ecov; y++) Ecov_x(y,i) = Ecov_mu + Ecov_re(y,i);
    }

    // add to nll if estimated (option in projection years to fix Ecov at last or average value)
    for(int y = 0; y < n_years_Ecov + n_years_proj_Ecov; y++){
      if(Ecov_use_re(y,i) == 1){
        nll += nll_Ecov(y,i);
      }
    }
  } // end loop over Ecovs
  SIMULATE if(simulate_state(3) == 1) if(sum(simulate_period) > 0) REPORT(Ecov_re);
  if(Ecov_model.sum() > 0) if(do_post_samp(3) == 1) ADREPORT(Ecov_re);

  // Environmental covariate observation model -------------------------------------
  //TODO: Ecov obs are not yet simulated in projection years!!!!!!!!
  Type nll_Ecov_obs = Type(0);
  Type nll_Ecov_obs_sig = Type(0); // Ecov obs sigma random effects (opt = 4)
  matrix<Type> Ecov_obs_sigma(n_years_Ecov, n_Ecov);
  for(int i = 0; i < n_Ecov; i++){
    for(int y = 0; y < n_years_Ecov; y++){
      if(Ecov_obs_sigma_opt(i) == 4){
        Type mu_logsigma = Ecov_obs_sigma_par(0,i);
        Type sd_logsigma = exp(Ecov_obs_sigma_par(1,i));
        nll_Ecov_obs_sig -= dnorm(Ecov_obs_logsigma_re(y,i), mu_logsigma, sd_logsigma, 1);
        SIMULATE if(simulate_data(2) == 1) if(simulate_period(0) == 1) {
          Ecov_obs_logsigma_re(y,i) = rnorm(mu_logsigma, sd_logsigma);
        }
        Ecov_obs_sigma(y,i) = exp(Ecov_obs_logsigma_re(y,i));
      } else{
        Ecov_obs_sigma(y,i) = exp(Ecov_obs_logsigma(y,i));
      }
      if(Ecov_use_obs(y,i) == 1){
        nll_Ecov_obs -= keep(keep_E(y,i)) * dnorm(obsvec(keep_E(y,i)), Ecov_x(y,i), Ecov_obs_sigma(y,i), 1);
        nll_Ecov_obs -= keep.cdf_lower(keep_E(y,i)) * log(squeeze(pnorm(obsvec(keep_E(y,i)), Ecov_x(y,i), Ecov_obs_sigma(y,i))));
        nll_Ecov_obs -= keep.cdf_upper(keep_E(y,i)) * log(1.0 - squeeze(pnorm(obsvec(keep_E(y,i)), Ecov_x(y,i), Ecov_obs_sigma(y,i))));
        SIMULATE if(simulate_data(2) ==1) if(simulate_period(0) == 1) {
          Ecov_obs(y,i) = rnorm(Ecov_x(y,i), Ecov_obs_sigma(y,i));
          obsvec(keep_E(y,i)) = Ecov_obs(y,i);
        }
      }
    }
  }
  nll += nll_Ecov_obs_sig;
  nll += nll_Ecov_obs;
  REPORT(nll_Ecov_obs_sig);
  SIMULATE if(simulate_data(2) ==1) if(simulate_period(0) == 1) {
    REPORT(Ecov_obs);
    REPORT(Ecov_obs_logsigma);
  }

 
  // Lag environmental covariates -------------------------------------
  // Then use Ecov_out(t) for processes in year t, instead of Ecov_x
  int n_effects = Ecov_beta.dim(0); // recruitment, mortality and any catchabilities, LAA, WAA, maturity
  array<Type> Ecov_out(n_years_model + n_years_proj, n_effects, n_Ecov); // Pop model uses Ecov_out(t) for processes in year t (Ecov_x shifted by lag and padded)
  Ecov_out.setZero(); // set Ecov_out = 0
  for(int i = 0; i < n_Ecov; i++){
    for(int t = 0; t < n_effects; t++){
      int ct = 0;
      for(int y = ind_Ecov_out_start(i,t); y < ind_Ecov_out_end(i,t) + 1 + n_years_proj; y++){
        Ecov_out(ct,t,i) = Ecov_x(y,i);
        ct++;
      }
    }
  }

  // Calculate ecov link model (b1*ecov + b2*ecov^2 + ...) --------------------
  // ecov_beta is now 4D array, dim = (2 + n_indices) x n_poly x n_ecov x n_ages
  int n_poly = Ecov_beta.dim(1); // now a 4D array dim: (n_effects,n_poly,n_Ecov,n_ages) is second dimension
  //vector<matrix<Type>> Ecov_lm(n_Ecov)(n_effects); // ecov linear model for each Ecov, dim = n_years_model + n_years_proj, n_ages
  // Ecov_lm.setZero();
  // Ecov_lm stores the linear models for each Ecov and where it is used. dim = n_Ecov, n_effects, n_years_model + n_years_proj, n_ages
  // n_effects dimension is: 0: recruitment, 1: M, 2-1+n_indices: which catchability it affects
  array<Type> Ecov_lm(n_Ecov, n_effects,n_years_model + n_years_proj, n_ages); 
  //vector<matrix<Type> > Ecov_lm(n_Ecov); // ecov linear model for each Ecov, dim = n_years_model + n_years_proj, n_ages
  for(int i = 0; i < n_Ecov; i++){
    for(int t = 0; t < n_effects; t++){
      vector<Type> thecol(n_years_model + n_years_proj);
      for(int y = 0; y < n_years_model + n_years_proj; y++) thecol(y) = Ecov_out(y,t,i);
      matrix<Type> X_poly(n_years_model + n_years_proj, n_poly);
      X_poly.setZero();
      if(n_poly == 1){ // n_poly = 1 if ecov effect is none or linear
        X_poly = thecol.matrix();
      } else { // n_poly > 1, get poly transformation for ith ecov
        X_poly = poly_trans(thecol, n_poly, n_years_model, n_years_proj);
      }
      for(int y = 0; y < n_years_model + n_years_proj; y++){
        for(int a = 0; a < n_ages; a++){
          for(int j = 0; j < n_poly; j++){
            Ecov_lm(i,t,y,a) += Ecov_beta(t,j,i,a) * X_poly(y,j); // poly transformation returns design matrix, don't need to take powers
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------
  // LAA random effects section
  Type nll_LAA = Type(0);
  vector<Type> sigma_LAA(n_LAA_dim);
  vector<Type> rho_LAA_a(n_LAA_dim); // correlation age
  vector<Type> rho_LAA_y(n_LAA_dim); // correlation year
  vector<Type> prho_LAA_a(n_LAA_dim); // partial correlation age 3D
  vector<Type> prho_LAA_y(n_LAA_dim); // partial correlation year 3D
  vector<Type> prho_LAA_c(n_LAA_dim); // partial correlation cohort 3D
  Type Sigma_LAA = Type(0);
  
  for(int j = 0; j < n_LAA_dim; j++) {
	  
	  if(LAA_re_model(j) > 1) { // only if random effects on LAA active
	  
	  	sigma_LAA(j) = exp(LAA_repars(j,0)); // first RE parameter
		rho_LAA_a(j) = rho_trans(LAA_repars(j,1));
		rho_LAA_y(j) = rho_trans(LAA_repars(j,2));  
		prho_LAA_a(j) = LAA_repars(j,1);
		prho_LAA_y(j) = LAA_repars(j,2);  
		prho_LAA_c(j) = LAA_repars(j,3);	
	  
		if((LAA_model == 1)|(LAA_model == 2)) { // vB-classic or Richards
			  
			  if((LAA_re_model(j) == 2) | (LAA_re_model(j) == 4)) { // Only for ii_y and Ar1_y.
				  
				// likelihood of LAA parameters deviations
				vector<Type> GWre0(n_years_model + n_years_proj);
				for(int y = 0; y < n_years_model + n_years_proj; y++) GWre0(y) = LAA_re(y,0,j); 
				Sigma_LAA = pow(pow(sigma_LAA(j),2) / (1-pow(rho_LAA_y(j),2)),0.5);
				nll_LAA += SCALE(AR1(rho_LAA_y(j)), Sigma_LAA)(GWre0);
				SIMULATE if(simulate_state(5) == 1) if(sum(simulate_period) > 0) {
					AR1(rho_LAA_y(j)).simulate(GWre0);
					for(int i = 0; i < GWre0.size(); i++) GWre0(i) = Sigma_LAA * GWre0(i);
					  for(int y = 0; y < n_years_model + n_years_proj; y++){
						if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
						  for(int a = 0; a < n_ages; a++) LAA_re(y,a,j) = GWre0(y);
						}
					  }
				}
				
			  } // LAA_re_model = 2 or 4
					  
			  if((LAA_re_model(j) == 3) | (LAA_re_model(j) == 5)) { // Only for ii_c and Ar1_c
				
				// likelihood of LAA parameters deviations				
				vector<Type> GWre0(n_years_model + n_years_proj + n_ages - 1);
				for(int i = 0; i < (n_ages - 1); i++) GWre0(i) = LAA_re(0,n_ages - i - 1,j); // for cohorts at y = 0 except a = 0
				for(int i = (n_ages - 1); i < (n_years_model + n_years_proj + n_ages - 1); i++) GWre0(i) = LAA_re(i-n_ages+1,0,j); // for cohorts y>=0 
				Sigma_LAA = pow(pow(sigma_LAA(j),2) / (1-pow(rho_LAA_y(j),2)),0.5);
				nll_LAA += SCALE(AR1(rho_LAA_y(j)), Sigma_LAA)(GWre0);
				SIMULATE if(simulate_state(5) == 1) if(sum(simulate_period) > 0) {
				  AR1(rho_LAA_y(j)).simulate(GWre0);
				  for(int i = 0; i < GWre0.size(); i++) GWre0(i) = Sigma_LAA * GWre0(i);
				  for(int y = 0; y < n_years_model + n_years_proj; y++){
					if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
					  for(int a = 0; a < n_ages; a++) {
						  if(y == 0) LAA_re(y,n_ages-1-a,j) = GWre0(a); 
						  else LAA_re(y,a,j) = GWre0(y-a+n_ages-1);
					  }
					}
				  }
				}
						  
			  } // LAA_re_model = 3 or 5
				  
		  } // If LAA_model = 1 or 2
  
  
		  if(LAA_model == 3) { // Nonparametric model
				  
				array<Type> LAAre0 = LAA_re.col(j);
				if((LAA_re_model(j) == 2) | (LAA_re_model(j) == 3)) { // iid and 2DAR1
				  // likelihood of LAA deviations
				  Sigma_LAA = pow(pow(sigma_LAA(j),2) / ((1-pow(rho_LAA_y(j),2))*(1-pow(rho_LAA_a(j),2))),0.5);
				  nll_LAA += SCALE(SEPARABLE(AR1(rho_LAA_a(j)),AR1(rho_LAA_y(j))), Sigma_LAA)(LAAre0); // must be array, not matrix!
				  SIMULATE if(simulate_state(5) == 1) if(sum(simulate_period) > 0) {
					//array<Type> LAAre_tmp = LAA_re;
					SEPARABLE(AR1(rho_LAA_a(j)),AR1(rho_LAA_y(j))).simulate(LAAre0);
					LAAre0 = Sigma_LAA * LAAre0;
					for(int y = 0; y < n_years_model + n_years_proj; y++){
					  if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
						for(int a = 0; a < n_ages; a++) LAA_re(y,a,j) = LAAre0(y,a);
					  }
					}
				  }
				} 
				
				if(LAA_re_model(j) == 4) { // 3D smoother. Author: Cheng et al. (10.1016/j.fishres.2023.106755)
				  // Define precision matrix for GMRF
				  int total_n = (n_years_model + n_years_proj)*n_ages;
				  Eigen::SparseMatrix<Type> Q_sparseG(total_n, total_n); // Precision matrix
				  // Construct precision matrix here
				  Q_sparseG = construct_Q(n_years_model+n_years_proj, n_ages, ay3D_IndexL, prho_LAA_y(j), prho_LAA_a(j), prho_LAA_c(j), sigma_LAA(j), Var3D_ParamL);
				  nll_LAA += GMRF(Q_sparseG)(LAAre0); 
				  SIMULATE if(simulate_state(5) == 1) if(sum(simulate_period) > 0) {
					vector<Type> LAAre_tmp(total_n); // should be a vector
					GMRF(Q_sparseG).simulate(LAAre_tmp);
					for(int y = 0; y < n_years_model + n_years_proj; y++){
					  if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
						for(int a = 0; a < n_ages; a++) LAA_re(y,a,j) = LAAre_tmp((n_years_model+n_years_proj)*a + y); // sort vector into WAA_re array
					  }
					}
				  }
				}		  
					  		  
		  } // If LAA_model = 3
	  
	  } // If LAA_re_model > 1
  
  } // Loop number of RE dim
  
	if(do_post_samp.sum()==0){
		ADREPORT(sigma_LAA);
		ADREPORT(rho_LAA_a);
		ADREPORT(rho_LAA_y);
		ADREPORT(prho_LAA_a);
		ADREPORT(prho_LAA_y);
		ADREPORT(prho_LAA_c);
	}
				
  // Report quantities:
  REPORT(nll_LAA);
  nll += nll_LAA; 
  REPORT(LAA_a);
  REPORT(LAA_re);
  REPORT(LAA_repars);
  if(do_post_samp(5) == 1) ADREPORT(LAA_re);


  // ---------------------------------------------------------------------
  // WAA random effects section
  Type nll_WAA = Type(0);
  vector<Type> sigma_WAA(n_WAA_dim);
  vector<Type> rho_WAA_a(n_WAA_dim); // correlation age
  vector<Type> rho_WAA_y(n_WAA_dim); // correlation year
  vector<Type> prho_WAA_a(n_WAA_dim); // partial correlation age 3D
  vector<Type> prho_WAA_y(n_WAA_dim); // partial correlation year 3D
  vector<Type> prho_WAA_c(n_WAA_dim); // partial correlation cohort 3D
  Type Sigma_WAA = Type(0);
  
  if(isW_parametric == 1) {
  
	  for(int j = 0; j < n_WAA_dim; j++) {
		  
		  if(WAA_re_model(j) > 1) { // only if random effects on LAA active
		  
			sigma_WAA(j) = exp(WAA_repars(j,0)); // first RE parameter
			rho_WAA_a(j) = rho_trans(WAA_repars(j,1));
			rho_WAA_y(j) = rho_trans(WAA_repars(j,2));  
			prho_WAA_a(j) = WAA_repars(j,1);
			prho_WAA_y(j) = WAA_repars(j,2);  
			prho_WAA_c(j) = WAA_repars(j,3);	
		  
			if(WAA_model == 1) { // Allometric
				  
				  if((WAA_re_model(j) == 2) | (WAA_re_model(j) == 3)) { // Only for ii_y and Ar1_y.
					  
					// likelihood of WAA parameters deviations
					vector<Type> Wre0(n_years_model + n_years_proj);
					for(int y = 0; y < n_years_model + n_years_proj; y++) Wre0(y) = WAA_re(y,0,j); 
					Sigma_WAA = pow(pow(sigma_WAA(j),2) / (1-pow(rho_WAA_y(j),2)),0.5);
					nll_WAA += SCALE(AR1(rho_WAA_y(j)), Sigma_WAA)(Wre0);
					SIMULATE if(simulate_state(6) == 1) if(sum(simulate_period) > 0) {
						AR1(rho_WAA_y(j)).simulate(Wre0);
						for(int i = 0; i < Wre0.size(); i++) Wre0(i) = Sigma_WAA * Wre0(i);
						  for(int y = 0; y < n_years_model + n_years_proj; y++){
							if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
							  for(int a = 0; a < n_ages; a++) WAA_re(y,a,j) = Wre0(y);
							}
						  }
					}
					
				  } // WAA_re_model = 2 or 3
					  
			  } // If WAA_model = 1 
	  
	  
			  if(WAA_model == 2) { // Nonparametric model
					  
					array<Type> WAAre0 = WAA_re.col(j);
					if((WAA_re_model(j) == 2) | (WAA_re_model(j) == 3)) { // iid and 2DAR1
					  // likelihood of LAA deviations
					  Sigma_WAA = pow(pow(sigma_WAA(j),2) / ((1-pow(rho_WAA_y(j),2))*(1-pow(rho_WAA_a(j),2))),0.5);
					  nll_WAA += SCALE(SEPARABLE(AR1(rho_WAA_a(j)),AR1(rho_WAA_y(j))), Sigma_WAA)(WAAre0); // must be array, not matrix!
					  SIMULATE if(simulate_state(6) == 1) if(sum(simulate_period) > 0) {
						//array<Type> WAAre_tmp = WAA_re;
						SEPARABLE(AR1(rho_WAA_a(j)),AR1(rho_WAA_y(j))).simulate(WAAre0);
						WAAre0 = Sigma_WAA * WAAre0;
						for(int y = 0; y < n_years_model + n_years_proj; y++){
						  if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
							for(int a = 0; a < n_ages; a++) WAA_re(y,a,j) = WAAre0(y,a);
						  }
						}
					  }
					} 
					
					if(WAA_re_model(j) == 4) { // 3D smoother. Author: Cheng et al. (10.1016/j.fishres.2023.106755)
					  // Define precision matrix for GMRF
					  int total_n = (n_years_model + n_years_proj)*n_ages;
					  Eigen::SparseMatrix<Type> Q_sparseW(total_n, total_n); // Precision matrix
					  // Construct precision matrix here
					  Q_sparseW = construct_Q(n_years_model+n_years_proj, n_ages, ay3D_IndexW, prho_WAA_y(j), prho_WAA_a(j), prho_WAA_c(j), sigma_WAA(j), Var3D_ParamW);
					  nll_WAA += GMRF(Q_sparseW)(WAAre0); 
					  SIMULATE if(simulate_state(6) == 1) if(sum(simulate_period) > 0) {
						vector<Type> WAAre_tmp(total_n); // should be a vector
						GMRF(Q_sparseW).simulate(WAAre_tmp);
						for(int y = 0; y < n_years_model + n_years_proj; y++){
						  if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
							for(int a = 0; a < n_ages; a++) WAA_re(y,a,j) = WAAre_tmp((n_years_model+n_years_proj)*a + y); // sort vector into WAA_re array
						  }
						}
					  }
					}		  
								  
			  } // If WAA_model = 2
		  
		  } // If WAA_re_model > 1
	  
	  } // Loop number of RE dim
  
  } // If W parametric
  
	if(do_post_samp.sum()==0){
		ADREPORT(sigma_WAA);
		ADREPORT(rho_WAA_a);
		ADREPORT(rho_WAA_y);
		ADREPORT(prho_WAA_a);
		ADREPORT(prho_WAA_y);
		ADREPORT(prho_WAA_c);
	}
				
  // Report quantities:
  REPORT(nll_WAA);
  nll += nll_WAA; 
  REPORT(WAA_a);
  REPORT(WAA_re);
  REPORT(WAA_repars);
  if(do_post_samp(6) == 1) ADREPORT(WAA_re);

  // --------------------------------------------------------------------------
  // Maturity random effects section:
  Type nll_mat = Type(0);
  
  // For all growth parameters:
  vector<Type> sigma_mat(n_mat_par); // first RE parameter
  vector<Type> rho_mat_y(n_mat_par);  // second RE parameter

  if(isMat_parametric == 1) {

	  for(int j = 0; j < n_mat_par; j++) {
	  
	  Type Sigma_mat = Type(0);

		  if((mat_re_model(j) == 2) | (mat_re_model(j) == 3)) { // Only for ii_y and Ar1_y. 
			
			sigma_mat(j) = exp(mat_repars(j,0));
			rho_mat_y(j) = rho_trans(mat_repars(j,1));  
			// likelihood of maturity parameters deviations
				vector<Type> matre0(n_years_model + n_years_proj);
				for(int y = 0; y < n_years_model + n_years_proj; y++) matre0(y) = mat_re(y,j); 
				Sigma_mat = pow(pow(sigma_mat(j),2) / (1-pow(rho_mat_y(j),2)),0.5);
				nll_mat += SCALE(AR1(rho_mat_y(j)), Sigma_mat)(matre0);
				SIMULATE if(simulate_state(7) == 1) if(sum(simulate_period) > 0) {
				  AR1(rho_mat_y(j)).simulate(matre0);
				  for(int i = 0; i < matre0.size(); i++) matre0(i) = Sigma_mat * matre0(i);
				  for(int y = 0; y < n_years_model + n_years_proj; y++){
					if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
					  for(int a = 0; a < n_ages; a++) mat_re(y,j) = matre0(y);
					}
				  }
				}

			if(do_post_samp.sum()==0){
			  ADREPORT(sigma_mat);
			  ADREPORT(rho_mat_y);
			}
			
		  }
	  
	  }
	  
  } // If isMat parametric
  
  REPORT(nll_mat);
  nll += nll_mat; 
  REPORT(mat_a);
  REPORT(mat_re);
  REPORT(mat_repars);
  if(do_post_samp(7) == 1) ADREPORT(mat_re);


  // ---------------------------------------------------------------------
  // LAA parameters per year,age
  array<Type> LAA_par(n_years_model + n_years_proj, n_ages, n_LAA_dim); // array for LAA parameters
  // Construct LAA parameters during model period:
  for(int j = 0; j < n_LAA_dim; j++) { 
	  for(int y = 0; y < n_years_model; y++) { 
			for(int a = 0; a < n_ages; a++) { 
				if((LAA_model == 1)|(LAA_model == 2)) LAA_par(y,a,j) = exp(LAA_a(j) + LAA_re(y,a,j)); // parametric options
				if(LAA_model == 3) LAA_par(y,a,j) = exp(LAA_a(a) + LAA_re(y,a,j)); // nonparametric 
			}
	  }
   }
  // add to LAA parameters in projection years
  if(do_proj == 1){  
	for(int j = 0; j < n_LAA_dim; j++) { 
	  if(proj_LAA_opt(j) == 2){ // average
		  matrix<Type> GW_toavg(n_toavg, n_ages);
		  for(int a = 0; a < n_ages; a++){
			for(int i = 0; i < n_toavg; i++){
			  GW_toavg(i,a) = LAA_par(avg_years_ind(i),a,j);
			}
		  }
		  vector<Type> GW_proj = GW_toavg.colwise().mean();
		  for(int y = n_years_model; y < n_years_model + n_years_proj; y++){
			for(int a = 0; a < n_ages; a++){
			  LAA_par(y,a,j) = GW_proj(a);
			}
		  }
	  } else { // proj_LAA_opt == 1
		for(int y = n_years_model; y < n_years_model + n_years_proj; y++) {
			for(int a = 0; a < n_ages; a++) { 		
				if((LAA_model == 1)|(LAA_model == 2)) LAA_par(y,a,j) = exp(LAA_a(j) + LAA_re(y,a,j)); // parametric options
				if(LAA_model == 3) LAA_par(y,a,j) = exp(LAA_a(a) + LAA_re(y,a,j)); // nonparametric 				
			}
		}
	  }
	}
  }
  // add ecov effect on LAA paramteres
  for(int j = 0; j < n_LAA_dim; j++) { 
	for(int i=0; i < n_Ecov; i++){
		if(Ecov_where(i,2+n_indices+j) == 1) {  
			for(int y = 0; y < n_years_model + n_years_proj; y++) {
				for(int a = 0; a < n_ages; a++) { 
					if((LAA_model == 1)|(LAA_model == 2)) LAA_par(y,a,j) *= exp(Ecov_lm(i,2+n_indices+j,y,0));
					if(LAA_model == 3) LAA_par(y,a,j) *= exp(Ecov_lm(i,2+n_indices+j,y,a)); // effect shared across ages
				}
			}
		}
	}
  }
  
  // Exp SD_par:
  vector<Type> SD_len(2); // always 2 parameters
  for(int i=0; i < 2; i++) {
	  SD_len(i) = exp(SD_par(i)); // parametric and semiparametric
  }
  REPORT(SD_len);

  // ---------------------------------------------------------------------
  // WAA parameters per year,age
  array<Type> WAA_par(n_years_model + n_years_proj, n_ages, n_WAA_dim); // array for LAA parameters
  // Construct WAA parameters during model period:
  for(int j = 0; j < n_WAA_dim; j++) { 
	  for(int y = 0; y < n_years_model; y++) { 
			for(int a = 0; a < n_ages; a++) { 
				if(WAA_model == 1) WAA_par(y,a,j) = exp(WAA_a(j) + WAA_re(y,a,j)); // parametric options
				if(WAA_model == 2) WAA_par(y,a,j) = exp(WAA_a(a) + WAA_re(y,a,j)); // nonparametric 
			}
	  }
   }
  // add to WAA parameters in projection years
  if(do_proj == 1){  
	for(int j = 0; j < n_WAA_dim; j++) {
	  if(proj_WAA_opt(j) == 2){ // average
		  matrix<Type> W_toavg(n_toavg, n_ages);
		  for(int a = 0; a < n_ages; a++){
			for(int i = 0; i < n_toavg; i++){
			  W_toavg(i,a) = WAA_par(avg_years_ind(i),a,j);
			}
		  }
		  vector<Type> W_proj = W_toavg.colwise().mean();
		  for(int y = n_years_model; y < n_years_model + n_years_proj; y++){
			for(int a = 0; a < n_ages; a++){
			  WAA_par(y,a,j) = W_proj(a);
			}
		  }
	  } else { // proj_WAA_opt == 1
		for(int y = n_years_model; y < n_years_model + n_years_proj; y++) {
			for(int a = 0; a < n_ages; a++) { 		
				if(WAA_model == 1) WAA_par(y,a,j) = exp(WAA_a(j) + WAA_re(y,a,j)); // parametric options
				if(WAA_model == 2) WAA_par(y,a,j) = exp(WAA_a(a) + WAA_re(y,a,j)); // nonparametric 				
			}
		}
	  }
	}
  }
  // add ecov effect on WAA paramteres
  for(int j = 0; j < n_WAA_dim; j++) { 
	for(int i=0; i < n_Ecov; i++){
		if(Ecov_where(i,2+n_indices+n_LAA_dim+j) == 1) {  
			for(int y = 0; y < n_years_model + n_years_proj; y++) {
				for(int a = 0; a < n_ages; a++) { 
					if(WAA_model == 1) WAA_par(y,a,j) *= exp(Ecov_lm(i,2+n_indices+n_LAA_dim+j,y,0));
					if(WAA_model == 2) WAA_par(y,a,j) *= exp(Ecov_lm(i,2+n_indices+n_LAA_dim+j,y,a)); // effect shared across ages
				}
			}
		}
	}
  }

  // ---------------------------------------------------------------------
  // Maturity per year
  matrix<Type> mat_par(n_years_model + n_years_proj,n_mat_par); // array for maturity parameters
  // Construct Maturity parameters during model period
  for(int j = 0; j < n_mat_par; j++) { 
	  for(int y = 0; y < n_years_model; y++) { 
			mat_par(y,j) = exp(mat_a(j) + mat_re(y,j)); 
	  }
   }
  // add to maturity parameters in projection years
  if(do_proj == 1){ 
	for(int j = 0; j < n_mat_par; j++) { 
	  if(proj_mat_opt(j) == 2){
		  vector<Type> mat_toavg(n_toavg);
		  for(int i = 0; i < n_toavg; i++){
			 mat_toavg(i) = mat_par(avg_years_ind(i),j);
		  }
		  Type mat_proj = mat_toavg.mean();
		  for(int y = n_years_model; y < n_years_model + n_years_proj; y++){
			  mat_par(y,j) = mat_proj;
		  }
	  } else { // proj_mat_opt == 1
		for(int y = n_years_model; y < n_years_model + n_years_proj; y++) {		
				mat_par(y,j) = exp(mat_a(j) + mat_re(y,j)); 
		}
	  }
	}
  }
  // add ecov effect on maturity paramteres
  for(int j = 0; j < n_mat_par; j++) { 
	for(int i=0; i < n_Ecov; i++){
		if(Ecov_where(i,2+n_indices+n_LAA_dim+n_WAA_dim+j) == 1) { 
			for(int y = 0; y < n_years_model + n_years_proj; y++) {
				mat_par(y,j) *= exp(Ecov_lm(i,2+n_indices+n_LAA_dim+n_WAA_dim+j,y,0));
			}
		}
	}
  }
  
  // --------------------------------------------------------------------------
  // Calculate mean-LAA, SDAA, and transition matrix, for all years (January 1st): 
  // This considers the random effects
  Type Slope = 0.0; // for SDAA calculation
  vector<Type> GW_fix_vector(n_LAA_par);
  // 1) Parametric approaches:
  if((LAA_model == 1)|(LAA_model == 2)) {
	GW_fix_vector = exp(LAA_a);
	jan1LAA = calculate_jan1LAA(LAA_par, n_ages, n_years_model, n_years_proj, lengths, LAA_model, age_L1, age_L1_ceil, GW_fix_vector); 
  }
  // 2) Nonparametric approach:
  if(LAA_model == 3) jan1LAA = LAA_par.col(0).matrix(); 

  // calculate SD and transition matrix:
  for(int y = 0; y < n_years_model + n_years_proj; y++) {	  
	  for(int a = 0; a < n_ages; a++) {
			// SD calculation: 
			if((LAA_model == 1)|(LAA_model == 2)) { // Parametric approaches
				// SD at age as function of mean LAA:
				if((a + 1.0) < age_L1) { // same as SD1
					SDAA(y,a) = SD_len(0); 
				} else { 
					if(a == (n_ages-1)) { // same as SDA
						SDAA(y,a) = SD_len(1);
					} else { // linear interpolation
						Slope = (SD_len(1) - SD_len(0))/(LAA_par(y,a,1)-LAA_par(y,a,2));
						SDAA(y,a) = SD_len(0) + Slope*(jan1LAA(y,a)-LAA_par(y,a,2));  
					}
				}
				// SD at age as function of age:
				// Slope = (SD_len(1) - SD_len(0))/(n_ages-1.0);
				// SDAA(y,a) = SD_len(0) + Slope*(a);  				
			}
			if(LAA_model == 3) { // Nonparametric
				// SD at age as function of mean LAA:
				Slope = (SD_len(1) - SD_len(0))/(jan1LAA(y,n_ages-1)-jan1LAA(y,0));
				SDAA(y,a) = SD_len(0) + Slope*(jan1LAA(y,a)-jan1LAA(y,0));  
				// SD at age as function of age:
				// Slope = (SD_len(1) - SD_len(0))/(n_ages-1.0);
				// SDAA(y,a) = SD_len(0) + Slope*(a);  				
			}  
		} // loop age
			
		jan1_phi_mat.col(y) = construct_phi_matrix(lengths, vector<Type>(jan1LAA.row(y)), vector<Type>(SDAA.row(y))); // fill the phi_matrix	
		
  } // loop year

  if(do_post_samp.sum()==0) if(isW_parametric == 1) ADREPORT(jan1LAA);// only when LAA is relevant

  // --------------------------------------------------------------------------
  // Calculate phi matrix jan1 constructed from fixed effects (ignores random effects)
  matrix<Type> fix_jan1LAA(n_years_model + n_years_proj, n_ages); // we do not really need the y index since it is constant over time
  array<Type> LAA_par_fix(n_years_model + n_years_proj, n_ages, n_LAA_dim);
  int n_yrs = n_years_model + n_years_proj;
  // Parametric approaches:
  if((LAA_model == 1)|(LAA_model == 2)) {   
	for(int j = 0; j < n_LAA_dim; j++) {
		for(int y = 0; y < n_years_model + n_years_proj; y++) for(int a = 0; a < n_ages; a++) LAA_par_fix(y,a,j) = exp(LAA_a(j));
	}
	fix_jan1LAA = calculate_jan1LAA(LAA_par_fix, n_ages, n_years_model, n_years_proj, lengths, LAA_model, age_L1, age_L1_ceil, GW_fix_vector); 
  }
  // Nonparametric
  if(LAA_model == 3) {
	for(int y = 0; y < n_years_model + n_years_proj; y++) fix_jan1LAA.row(y) = exp(LAA_a); // nonparametric approach
  }
  // Fixed phi matrix jan-1: This should be used to calculate reference points?
  array<Type> fix_phi_mat = fryr_phi_matrix(fix_jan1LAA, n_yrs, n_years_model, SD_len, LAA_par_fix, lengths, fracyr_catch, LAA_model, age_L1);; // 
  REPORT(fix_phi_mat);

  // --------------------------------------------------------------------------
  // Construct phi matrices:
  // Do it here just once (more efficient?)
  array<Type> out_phi_mat(n_lengths, n_ages, n_years_model + n_years_proj); // used later
  array<Type> ssb_phi_mat(n_lengths, n_ages, n_years_model + n_years_proj); // used later for SSB
  array<Type> catch_phi_mat(n_lengths, n_ages, n_years_model + n_years_proj); // used later for catch
  // January 1st:
  phi_matrix(waa_pointer_jan1-1) = jan1_phi_mat; // calculated in previous section
  // SSB:
  phi_matrix(waa_pointer_ssb-1) = fryr_phi_matrix(jan1LAA, n_yrs, n_years_model, SD_len, LAA_par, lengths, fracyr_SSB, LAA_model, age_L1);
  // Total catch:
  phi_matrix(waa_pointer_totcatch-1) = fryr_phi_matrix(jan1LAA, n_yrs, n_years_model, SD_len, LAA_par, lengths, fracyr_catch, LAA_model, age_L1);
  // For fleets:
  for(int f = 0; f < n_fleets; f++) {
	 phi_matrix(waa_pointer_fleets(f)-1) = phi_matrix(waa_pointer_totcatch-1); // use same as total catch, fraction = 0.5
  }
  // For indices:
  for(int i = 0; i < n_indices; i++) {
	 phi_matrix(waa_pointer_indices(i)-1) = fryr_phi_matrix(jan1LAA, n_yrs, n_years_model, SD_len, LAA_par, lengths, vector<Type>(fracyr_indices.col(i)), LAA_model, age_L1);
  }
  
  // --------------------------------------------------------------------------
  // Weight at age calculations:
  // Exclusively for Allometric and used for SSB calculations
  matrix<Type> wt_at_age(n_years_model + n_years_proj, n_ages);
  matrix<Type> wt_at_len(n_years_model + n_years_proj, n_lengths);
  wt_at_len.setZero();
  wt_at_age.setZero();
  Type sum_wt_ssb = 0;
  Type lenmid = (lengths(1) - lengths(0))*0.5;

  // weight (Allometric approach):
  if((isW_parametric == 1) & (WAA_model == 1)) { 
	// First calculate weight at length:
	for(int y = 0; y < n_years_model + n_years_proj; y++)  {
		for(int l = 0; l < n_lengths; l++) {
			wt_at_len(y,l) = WAA_par(y,0,0)*pow((lengths(l)+lenmid), WAA_par(y,0,1)); // age = 0
		}
	}  
	// Then calculate weight at age:
	out_phi_mat = phi_matrix(waa_pointer_ssb-1);
	for(int y = 0; y < n_years_model + n_years_proj; y++)  {
		for(int a = 0; a < n_ages; a++) {
			sum_wt_ssb = 0;
			for(int l = 0; l < n_lengths; l++) {
				sum_wt_ssb += out_phi_mat(l,a,y)*wt_at_len(y,l);
			}
			wt_at_age(y,a) = sum_wt_ssb;
		}
	}
  }
  REPORT(wt_at_len);
  REPORT(wt_at_age);
  
  // --------------------------------------------------------------------------
  // Calculate maturity at age regardless the chosen method: 
  // Used in SSB calculations  
  matrix<Type> mat_at_age(n_years_model + n_years_proj, n_ages);
  matrix<Type> mat_at_len(n_years_model + n_years_proj, n_lengths);
  mat_at_len.setZero();
  Type sum_mat = 0;
  
  if(isMat_parametric == 0) { 
	mat_at_age = mature; // same as data input
  }
  if((isMat_parametric == 1) & (mat_model == 1)) { // age-logistic parametric approach
	for(int y = 0; y < n_years_model + n_years_proj; y++)  for(int a = 0; a < n_ages; a++) mat_at_age(y,a) = 1/(1+exp(-mat_par(y,0)*((a+1) - mat_par(y,1))));
  }
  if((isMat_parametric == 1) & (mat_model == 2)) { // len-logistic parametric approach
	// First calculate maturity at length:
	for(int y = 0; y < n_years_model + n_years_proj; y++)  {
		for(int l = 0; l < n_lengths; l++) {
			mat_at_len(y,l) = 1/(1+exp(-mat_par(y,0)*((lengths(l)+lenmid) - mat_par(y,1))));
		}
	}
	// Then calculate maturity at age:
	out_phi_mat = phi_matrix(waa_pointer_ssb-1);
	for(int y = 0; y < n_years_model + n_years_proj; y++)  {
		for(int a = 0; a < n_ages; a++) {
			sum_mat = 0;
			for(int l = 0; l < n_lengths; l++) {
				sum_mat += out_phi_mat(l,a,y)*mat_at_len(y,l);
			}
			mat_at_age(y,a) = sum_mat;
		}
	}
  } 
  REPORT(mat_at_len);

  // --------------------------------------------------------------------------
  // Weight at age calculations:
	Type sum_wt = 0;
	Type sum_wt_fleet = 0;
	Type sum_wt_index = 0;
	vector<Type> t_pred_waa(n_ages);
	vector<Type> t_obs_waa(n_ages);
	vector<Type> t_cv_waa(n_ages);
	array<Type> waa_proj(waa.dim(0), n_years_proj, n_ages);
	array<Type> nll_waa(waa.dim(0), n_years_model);
	nll_waa.setZero();
  if(isW_parametric == 0) { // empirical weight-at-age
  	// Replace pred_waa by waa to be used later:
	for(int y = 0; y < n_years_model + n_years_proj; y++) {
		for(int a = 0; a < n_ages; a++) { 
			pred_waa(waa_pointer_jan1 - 1,y,a) = waa(waa_pointer_jan1 - 1,y,a);
			pred_waa(waa_pointer_ssb - 1,y,a) = waa(waa_pointer_ssb - 1,y,a);
			pred_waa(waa_pointer_totcatch - 1,y,a) = waa(waa_pointer_totcatch - 1,y,a);
			for(int f = 0; f < n_fleets; f++) {
				pred_waa(waa_pointer_fleets(f) - 1,y,a) = waa(waa_pointer_fleets(f)-1,y,a);
			}
			for(int i = 0; i < n_indices; i++) {
				pred_waa(waa_pointer_indices(i) - 1,y,a) = waa(waa_pointer_indices(i)-1,y,a);
			}
		}
	}
  } else {
		if(WAA_model == 1) { // Allometric
			  // For Jan-1
			  out_phi_mat = phi_matrix(waa_pointer_jan1-1);
			  for(int y = 0; y < n_years_model + n_years_proj; y++) { // 
				for(int a = 0; a < n_ages; a++) { 
					sum_wt = 0;
					for(int l = 0; l < n_lengths; l++) {
						sum_wt += out_phi_mat(l,a,y)*wt_at_len(y,l);
					}
					pred_waa(waa_pointer_jan1 - 1,y,a) = sum_wt; // jan-1st
				}
			  }
				// For SSB
			  for(int y = 0; y < n_years_model + n_years_proj; y++) { // 
				for(int a = 0; a < n_ages; a++) { 
					pred_waa(waa_pointer_ssb - 1,y,a) = wt_at_age(y,a); // SSB
				}
			  }
			  
				// For fleets
				for(int f = 0; f < n_fleets; f++) {
					out_phi_mat = phi_matrix(waa_pointer_fleets(f)-1);
					for(int y = 0; y < n_years_model + n_years_proj; y++) { // 
						for(int a = 0; a < n_ages; a++) { 
							sum_wt_fleet = 0;
							for(int l = 0; l < n_lengths; l++) {
								sum_wt_fleet += out_phi_mat(l,a,y)*wt_at_len(y,l);
							}
							pred_waa(waa_pointer_fleets(f)-1,y,a) = sum_wt_fleet; 
							pred_waa(waa_pointer_totcatch-1,y,a) = sum_wt_fleet;
							t_pred_waa(a) = sum_wt_fleet; // save it as vector predictions
						}
						if(y < n_years_model) if(use_catch_waa(y,f) == 1) { 
							t_obs_waa = obsvec.segment(keep_Cwaa(f,y,0), keep_Cwaa(f,y,1));	
							for(int a = 0; a < n_ages; a++) t_cv_waa(a) = waa_cv(waa_pointer_fleets(f) - 1,y,a);
							nll_waa(waa_pointer_fleets(f) - 1,y) -= get_waa_ll(t_obs_waa, t_pred_waa, t_cv_waa, bias_correct_oe);
						}
						SIMULATE if(simulate_data(0) == 1) if(use_catch_waa(y,f) == 1){
							if((simulate_period(0) == 1) & (y < n_years_model)) //model years
							{
								vector<Type> tf_waa_obs = sim_waa(t_pred_waa, t_cv_waa, bias_correct_pe);
								obsvec.segment(keep_Cwaa(f,y,0),keep_Cwaa(f,y,1)) = tf_waa_obs;
								for(int a = 0; a < n_ages; a++) waa(waa_pointer_fleets(f) - 1,y,a) = tf_waa_obs(a);
							}
							if((simulate_period(1) == 1) & (y > n_years_model - 1)) //projection years
							{
								for(int a = 0; a < n_ages; a++) t_cv_waa(a) = waa_cv(waa_pointer_fleets(f) - 1,n_years_model-1,a);// use last year CV
								vector<Type> tf_waa_obs = sim_waa(t_pred_waa, t_cv_waa, bias_correct_pe);
								for(int a = 0; a < n_ages; a++) waa_proj(waa_pointer_fleets(f) - 1,y,a) = tf_waa_obs(a);
							}
						}
					}
				}
				
				// For indices
				for(int i = 0; i < n_indices; i++) {
					out_phi_mat = phi_matrix(waa_pointer_indices(i)-1);
					for(int y = 0; y < n_years_model + n_years_proj; y++) { // 
						for(int a = 0; a < n_ages; a++) { 
							sum_wt_index = 0;
							for(int l = 0; l < n_lengths; l++) {
								sum_wt_index += out_phi_mat(l,a,y)*wt_at_len(y,l);
							}
							pred_waa(waa_pointer_indices(i)-1,y,a) = sum_wt_index; // for indices
							t_pred_waa(a) = sum_wt_index; // save it as vector predictions							
						}
						if(y < n_years_model) if(use_index_waa(y,i) == 1) { 
							t_obs_waa = obsvec.segment(keep_Iwaa(i,y,0), keep_Iwaa(i,y,1));	
							for(int a = 0; a < n_ages; a++) t_cv_waa(a) = waa_cv(waa_pointer_indices(i) - 1,y,a);
							nll_waa(waa_pointer_indices(i) - 1,y) -= get_waa_ll(t_obs_waa, t_pred_waa, t_cv_waa, bias_correct_oe);
						}
						SIMULATE if(simulate_data(1) == 1) if(use_index_waa(y,i) == 1){
							if((simulate_period(0) == 1) & (y < n_years_model)) //model years
							{
								vector<Type> tf_waa_obs = sim_waa(t_pred_waa, t_cv_waa, bias_correct_pe);
								obsvec.segment(keep_Iwaa(i,y,0),keep_Iwaa(i,y,1)) = tf_waa_obs;
								for(int a = 0; a < n_ages; a++) waa(waa_pointer_indices(i) - 1,y,a) = tf_waa_obs(a);
							}
							if((simulate_period(1) == 1) & (y > n_years_model - 1)) //projection years
							{
								for(int a = 0; a < n_ages; a++) t_cv_waa(a) = waa_cv(waa_pointer_indices(i) - 1,n_years_model-1,a);// use last year CV
								vector<Type> tf_waa_obs = sim_waa(t_pred_waa, t_cv_waa, bias_correct_pe);
								for(int a = 0; a < n_ages; a++) waa_proj(waa_pointer_indices(i) - 1,y,a) = tf_waa_obs(a);
							}
						}							
					}
				}
				
			  nll += nll_waa.sum();	
		} 
		
		if(WAA_model == 2) {// nonparametric approach
			vector<Type> fracyr_WAA(n_ages);
			matrix<Type> WAA_par_matrix(n_years_model + n_years_proj, n_ages);
			for(int y = 0; y < n_years_model + n_years_proj; y++) for(int a = 0; a < n_ages; a++) WAA_par_matrix(y,a) = WAA_par(y,a,0);
			for(int y = 0; y < n_years_model + n_years_proj; y++) {
				int yuse = y;
				int y_1 = y + 1;
				if(y > n_years_model - 1) yuse = n_years_model -1; //some things only go up to n_years_model-1
				if(y == (n_years_model + n_years_proj - 1)) y_1 = y;
				
				// For Jan-1
				for(int a = 0; a < n_ages; a++) { 
					pred_waa(waa_pointer_jan1 - 1,y,a) = WAA_par(y,a,0); // jan-1st
				}
				// For SSB
				fracyr_WAA = get_fracyr_WAA(vector<Type>(WAA_par_matrix.row(y)), vector<Type>(WAA_par_matrix.row(y_1)), fracyr_SSB(yuse));
				for(int a = 0; a < n_ages; a++) { 
					pred_waa(waa_pointer_ssb - 1,y,a) = fracyr_WAA(a); // SSB
				}	
				// For fleets
				for(int f = 0; f < n_fleets; f++) {
					fracyr_WAA = get_fracyr_WAA(vector<Type>(WAA_par_matrix.row(y)), vector<Type>(WAA_par_matrix.row(y_1)), fracyr_catch(yuse));
					for(int a = 0; a < n_ages; a++) { 
						pred_waa(waa_pointer_fleets(f)-1,y,a) = fracyr_WAA(a); 
						pred_waa(waa_pointer_totcatch-1,y,a) = fracyr_WAA(a); // for total catch, it is using the last fracyr_fleets
						t_pred_waa(a) = fracyr_WAA(a); // save it as vector predictions
					}						
					if(y < n_years_model) if(use_catch_waa(y,f) == 1) { 
						t_obs_waa = obsvec.segment(keep_Cwaa(f,y,0), keep_Cwaa(f,y,1));	
						for(int a = 0; a < n_ages; a++) t_cv_waa(a) = waa_cv(waa_pointer_fleets(f) - 1,y,a);
						nll_waa(waa_pointer_fleets(f) - 1,y) -= get_waa_ll(t_obs_waa, t_pred_waa, t_cv_waa, bias_correct_oe);
					}
					SIMULATE if(simulate_data(0) == 1) if(use_catch_waa(y,f) == 1){
						if((simulate_period(0) == 1) & (y < n_years_model)) //model years
						{
							vector<Type> tf_waa_obs = sim_waa(t_pred_waa, t_cv_waa, bias_correct_pe);
							obsvec.segment(keep_Cwaa(f,y,0),keep_Cwaa(f,y,1)) = tf_waa_obs;
							for(int a = 0; a < n_ages; a++) waa(waa_pointer_fleets(f) - 1,y,a) = tf_waa_obs(a);
						}
						if((simulate_period(1) == 1) & (y > n_years_model - 1)) //projection years
						{
							for(int a = 0; a < n_ages; a++) t_cv_waa(a) = waa_cv(waa_pointer_fleets(f) - 1,n_years_model-1,a);// use last year CV
							vector<Type> tf_waa_obs = sim_waa(t_pred_waa, t_cv_waa, bias_correct_pe);
							for(int a = 0; a < n_ages; a++) waa_proj(waa_pointer_fleets(f) - 1,y,a) = tf_waa_obs(a);
						}
					}	
				}
				
				// For indices
				for(int i = 0; i < n_indices; i++) {
					fracyr_WAA = get_fracyr_WAA(vector<Type>(WAA_par_matrix.row(y)), vector<Type>(WAA_par_matrix.row(y_1)), fracyr_indices(yuse,i));
					for(int a = 0; a < n_ages; a++) { 
						pred_waa(waa_pointer_indices(i)-1,y,a) = fracyr_WAA(a); // for indices	
						t_pred_waa(a) = fracyr_WAA(a); // save it as vector predictions							
					}
					if(y < n_years_model) if(use_index_waa(y,i) == 1) { 
						t_obs_waa = obsvec.segment(keep_Iwaa(i,y,0), keep_Iwaa(i,y,1));	
						for(int a = 0; a < n_ages; a++) t_cv_waa(a) = waa_cv(waa_pointer_indices(i) - 1,y,a);
						nll_waa(waa_pointer_indices(i) - 1,y) -= get_waa_ll(t_obs_waa, t_pred_waa, t_cv_waa, bias_correct_oe);
					}
					SIMULATE if(simulate_data(1) == 1) if(use_index_waa(y,i) == 1){
						if((simulate_period(0) == 1) & (y < n_years_model)) //model years
						{
							vector<Type> tf_waa_obs = sim_waa(t_pred_waa, t_cv_waa, bias_correct_pe);
							obsvec.segment(keep_Iwaa(i,y,0),keep_Iwaa(i,y,1)) = tf_waa_obs;
							for(int a = 0; a < n_ages; a++) waa(waa_pointer_indices(i) - 1,y,a) = tf_waa_obs(a);
						}
						if((simulate_period(1) == 1) & (y > n_years_model - 1)) //projection years
						{
							for(int a = 0; a < n_ages; a++) t_cv_waa(a) = waa_cv(waa_pointer_indices(i) - 1,n_years_model-1,a);// use last year CV
							vector<Type> tf_waa_obs = sim_waa(t_pred_waa, t_cv_waa, bias_correct_pe);
							for(int a = 0; a < n_ages; a++) waa_proj(waa_pointer_indices(i) - 1,y,a) = tf_waa_obs(a);
						}
					}							
				}				

			}
			nll += nll_waa.sum();	
		} // if nonparametric approach 		
  } // else isW_parametric
  REPORT(pred_waa);	
  REPORT(nll_waa);
  if(do_post_samp.sum()==0) if(isW_parametric == 1) ADREPORT(pred_waa); // If smoothing the WAA matrix get SEs
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(0) == 1) REPORT(waa); // fisheries
  SIMULATE if(simulate_data(1) == 1) if(simulate_period(0) == 1) REPORT(waa); // indices
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(1) == 1) REPORT(waa_proj); // fisheries
  SIMULATE if(simulate_data(1) == 1) if(simulate_period(1) == 1) REPORT(waa_proj); // indices


  // --------------------------------------------------------------------------
  // Calculate mortality (M, F, then Z)
  // Natural mortality process model
  Type nll_M = Type(0);
  if(M_re_model > 1) // random effects on M, M_re = 2D AR1 deviations on M(year,age), dim = n_years x n_M_a
  {
    Type sigma_M = exp(M_repars(0));
    Type rho_M_a = rho_trans(M_repars(1));
    Type rho_M_y = rho_trans(M_repars(2));
    Type Sigma_M;
    // likelihood of M deviations, M_re
    if((M_re_model == 2) | (M_re_model == 5)){ //2D AR1: age, year
      Sigma_M = pow(pow(sigma_M,2) / ((1-pow(rho_M_y,2))*(1-pow(rho_M_a,2))),0.5);
      nll_M += SCALE(SEPARABLE(AR1(rho_M_a),AR1(rho_M_y)), Sigma_M)(M_re); // must be array, not matrix!
      SIMULATE if(simulate_state(1) == 1) if(sum(simulate_period) > 0) {
        array<Type> Mre_tmp = M_re;
        SEPARABLE(AR1(rho_M_a),AR1(rho_M_y)).simulate(Mre_tmp);
        Mre_tmp = Sigma_M * Mre_tmp;
        for(int y = 0; y < n_years_model + n_years_proj; y++){
          if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
            for(int a = 0; a < n_M_a; a++) M_re(y,a) = Mre_tmp(y,a);
          }
        }
      }
    } else {
      if(M_re_model == 3){ // 1D ar1_a
        vector<Type> Mre0 = M_re.matrix().row(0);
        Sigma_M = pow(pow(sigma_M,2) / (1-pow(rho_M_a,2)),0.5);
        nll_M += SCALE(AR1(rho_M_a), Sigma_M)(Mre0);
        SIMULATE if(simulate_state(1) == 1) if(sum(simulate_period) > 0) {
          AR1(rho_M_a).simulate(Mre0);
          for(int i = 0; i < Mre0.size(); i++) Mre0(i) = Sigma_M * Mre0(i);
          for(int y = 0; y < n_years_model + n_years_proj; y++){
            for(int i = 0; i < Mre0.size(); i++){
              M_re(y,i) = Mre0(i);
            }
          }
        }          
      } else { // M_re_model = 4, 1D ar1_y
        vector<Type> Mre0 = M_re.matrix().col(0);
        Sigma_M = pow(pow(sigma_M,2) / (1-pow(rho_M_y,2)),0.5);
        nll_M += SCALE(AR1(rho_M_y), Sigma_M)(Mre0);
        SIMULATE if(simulate_state(1) == 1) if(sum(simulate_period) > 0) {
          AR1(rho_M_y).simulate(Mre0);
          for(int i = 0; i < Mre0.size(); i++) Mre0(i) = Sigma_M * Mre0(i);
          for(int y = 0; y < n_years_model + n_years_proj; y++){
            if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
              for(int a = 0; a < n_M_a; a++) M_re(y,a) = Mre0(y); //all ages mapped to the same annual RE
            }
          }
        }
      }
    }
    if(do_post_samp.sum()==0){
      ADREPORT(sigma_M);
      ADREPORT(rho_M_a);
      ADREPORT(rho_M_y);
    }
  }
  REPORT(nll_M);
  nll += nll_M;
  REPORT(M_re); //even if M_re not simulated.
  if(do_post_samp(1) == 1) ADREPORT(M_re);
  REPORT(M_a);
  REPORT(M_repars);

  // Construct mortality-at-age (MAA)
  matrix<Type> MAA(n_years_model + n_years_proj,n_ages);
  if(M_model == 2){ // age-specific M
    for(int a = 0; a < n_ages; a++) for(int y = 0; y < n_years_model; y++) MAA(y,a) = exp(M_a(a) + M_re(y,a));   
  } else {
    if(M_model == 1){ // constant M
      for(int a = 0; a < n_ages; a++) for(int y = 0; y < n_years_model; y++) MAA(y,a) = exp(M_a(0) + M_re(y,a));
    } else { // M_model = 3, M is allometric function of weight
      for(int a = 0; a < n_ages; a++) for(int y = 0; y < n_years_model; y++) MAA(y,a) = exp(M_a(0) + M_re(y,a) - exp(log_b) * log(pred_waa(waa_pointer_jan1-1,y,a)));
    }
  }
  // add to MAA in projection years
  if(do_proj == 1){ 
    if(proj_M_opt == 2){ // use average MAA over avg.yrs 
      matrix<Type> MAA_toavg(n_toavg,n_ages);
      for(int a = 0; a < n_ages; a++){
        for(int i = 0; i < n_toavg; i++){
          MAA_toavg(i,a) = MAA(avg_years_ind(i),a);
        }
      }
      vector<Type> MAA_proj = MAA_toavg.colwise().mean();
      for(int y = n_years_model; y < n_years_model + n_years_proj; y++){
        MAA.row(y) = MAA_proj;
      }
    } else { // proj_M_opt == 1, use M_re and/or ecov_lm in projection years
        if(M_model == 2){ // age-specific M
          for(int a = 0; a < n_ages; a++) for(int y = n_years_model; y < n_years_model + n_years_proj; y++) MAA(y,a) = exp(M_a(a) + M_re(y,a));   
        } else {
          if(M_model == 1){ // constant M
            for(int a = 0; a < n_ages; a++) for(int y = n_years_model; y < n_years_model + n_years_proj; y++) MAA(y,a) = exp(M_a(0) + M_re(y,a));
          } else { // M_model = 3, M is allometric function of weight
            for(int a = 0; a < n_ages; a++) for(int y = n_years_model; y < n_years_model + n_years_proj; y++) MAA(y,a) = exp(M_a(0) + M_re(y,a) - exp(log_b) * log(pred_waa(waa_pointer_jan1-1,y,a)));
          }
        }
      }
  }
  // add ecov effect on M (by year, shared across ages)
  for(int i=0; i < n_Ecov; i++){
    if(Ecov_where(i,1) == 1) { //not sure why Ecov_how is needed for M. if(Ecov_how(i) == 1){ // if ecov i affects M
      for(int a = 0; a < n_ages; a++){
        for(int y = 0; y < n_years_model + n_years_proj; y++) MAA(y,a) *= exp(Ecov_lm(i,1,y,a));
      }
    }
  }
  // prior on M(WAA) coefficient
  if(use_b_prior == 1)
  {
    Type mu = log(0.305);
    if(bias_correct_pe == 1) mu -= 0.5*exp(2*log(0.08));
    Type lprior_b = dnorm(log_b, mu, Type(0.08), 1);
    SIMULATE
    {
      if(simulate_state(1) == 1) if(sum(simulate_period) > 0) log_b = rnorm(mu, Type(0.08));
      REPORT(log_b);
    }
    REPORT(lprior_b);
    nll -= lprior_b;
  }

  // Survey catchability models
  matrix<Type> nll_q(n_years_model+n_years_proj,n_indices);
  nll_q.setZero();
  vector<Type> sigma_q(n_indices);
  sigma_q.setZero();
  vector<Type> rho_q(n_indices);
  rho_q.setZero();
  vector<Type> nll_q_prior(n_indices);
  nll_q_prior.setZero();
  matrix<Type> logit_q_mat(n_years_model+n_years_proj, n_indices);
  logit_q_mat.setZero();
  for(int i = 0; i < n_indices; i++) {
    
    //use prior for q? q_prior_re are random effects with mean logit_q (fixed) and sd = logit_q_prior_sigma.
    if(use_q_prior(i) == 1){ 
      nll_q_prior(i) -= dnorm(q_prior_re(i), logit_q(i), logit_q_prior_sigma(i), 1);
      SIMULATE if(simulate_state(4) == 1) if(sum(simulate_period) > 0){
        q_prior_re(i) = rnorm(logit_q(i), logit_q_prior_sigma(i));
      } 
      for(int y = 0; y < n_years_model + n_years_proj; y++) logit_q_mat(y,i) += q_prior_re(i);
    }
    else for(int y = 0; y < n_years_model + n_years_proj; y++) logit_q_mat(y,i) += logit_q(i);
    
    if(use_q_re(i) > 0) // random effects on q, q_re = AR1 deviations on (year,age), dim = n_years x n_M_a
    {
      sigma_q(i) = exp(q_repars(i,0)); // conditional sd
      rho_q(i) = rho_trans(q_repars(i,1)); // autocorrelation

      nll_q(0,i) -= dnorm(q_re(0,i), Type(0), sigma_q(i)*exp(-0.5 * log(1 - pow(rho_q(i),Type(2)))), 1);
      SIMULATE if((simulate_state(4) == 1) & (simulate_period(0) == 1)) {
        q_re(0,i) = rnorm(Type(0), sigma_q(i)*exp(-0.5 * log(1 - pow(rho_q(i),Type(2)))));
      }
      logit_q_mat(0,i) += q_re(0,i); //add in q random effects.
      for(int y = 1; y < n_years_model + n_years_proj; y++)
      {
        nll_q(y,i) -= dnorm(q_re(y,i), rho_q(i) * q_re(y-1,i), sigma_q(i), 1);
        SIMULATE if(simulate_state(4) == 1) {
          if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))) {
            q_re(y,i) = rnorm(rho_q(i) * q_re(y-1,i), sigma_q(i));
          }
        }
        logit_q_mat(y,i) += q_re(y,i); //add in q random effects.
      }
    }
  }

  nll += nll_q.sum() + nll_q_prior.sum();
  int Ecov_effects_on_q = 0;
  for(int y = 0; y < n_years_model + n_years_proj; y++) {
    for(int ind = 0; ind < n_indices; ind++) {
      for(int i=0; i < n_Ecov; i++){
        if(Ecov_where(i,2+ind) == 1){ // if ecov i affects q and which index
          logit_q_mat(y,ind) += Ecov_lm(i,2+ind,y,0);
          Ecov_effects_on_q++;
        }
      }
      q(y,ind) = q_lower(ind) + (q_upper(ind) - q_lower(ind))/(1 + exp(-logit_q_mat(y,ind)));
    }
  }
  
  // Transform selex-at-len to selex-at-age (only to save info in selAA), use JAN-1 phi matrix:
  // Selex-at-age is mandatory (for FAA calculation) even if selex-at-length is the main selectivity function
  // When transforming from selex-at-age to selex-at-length only put 1 for all length bins (placeholder)
  for(int b = 0; b < n_selblocks; b++) {
	  if(selblock_models(b) < 6) { // for selex-at-age models
		// selAA same as selAL (selex-at-age)
		selAA(b) = selAL(b);  
		// selLL = 1 (not important, just a placeholder)
		matrix<Type> tmpLL(n_years_model, n_lengths);
		for(int y = 0; y < n_years_model; y++) for(int l = 0; l < n_lengths; l++) tmpLL(y,l) = 1.0;
		selLL(b) = tmpLL;
	  } else { // for selex-at-length models
	  	// selLL same as selAL (selex-at-len)
		selLL(b) = selAL(b);
		// selLL = 1 (not important, just a placeholder)
		matrix<Type> tmpAA(n_years_model, n_ages);
		for(int y = 0; y < n_years_model; y++) for(int a = 0; a < n_ages; a++) tmpAA(y,a) = 1.0;
		selAA(b) = tmpAA;
	  }
  }
  
  // Construct survey catchability-at-age (QAA)
  // THIS SECTION IS NO LONGER REQUIRED 
  // for(int i = 0; i < n_indices; i++)
  // {
	// out_phi_mat = phi_matrix(waa_pointer_indices(i)-1);
  // // add ecov effect on M (by year, shared across ages)
    // for(int y = 0; y < n_years_model; y++)
    // {
	  // // It transforms selAL to selAA if required using phi_matrix at fracyr.
	  // matrix<Type> this_selAL = selAL(selblock_pointer_indices(y,i)-1);
	  // int this_sel_model = selblock_models(selblock_pointer_indices(y,i)-1);
	  // temp_selAA = get_selAA_from_selAL(this_selAL, y, this_sel_model, out_phi_mat);
      // for(int a = 0; a < n_ages; a++) {
		  // QAA(y,i,a) = q(y,i) * temp_selAA(a);
	  // }
    // }
    // //just use last years selectivity for now
    // if(do_proj == 1) {
		// for(int y = n_years_model; y < n_years_model + n_years_proj; y++) { 
			// for(int a = 0; a < n_ages; a++) {
				// QAA(y,i,a) = q(y,i) * temp_selAA(a); // using temp_selAA from last year (n_years_model)
			// }
		// }
	// }
  // }
  REPORT(logit_q_mat);
  if(use_q_re.sum()>0 || Ecov_effects_on_q>0) if(do_post_samp.sum()< 1) ADREPORT(logit_q_mat);
  REPORT(sigma_q);
  REPORT(rho_q);
  REPORT(nll_q);
  REPORT(nll_q_prior);
  REPORT(q_prior_re); //even if q_prior_re not simulated
  REPORT(q_re);
  if(do_post_samp(4)==1) ADREPORT(q_re); //even if q_re not simulated.
  REPORT(q);
  //REPORT(QAA);

  // Construct fishing mortality-at-age (FAA)
  FAA_tot.setZero();
  for(int f = 0; f < n_fleets; f++)
  {
    log_F(0,f) = log_F1(f);
    F(0,f) = exp(log_F(0,f));
	  // It transforms selAL to selAA if required using phi_matrix at fracyr.
	  catch_phi_mat = phi_matrix(waa_pointer_totcatch-1);
	  matrix<Type> this_selAL = selAL(selblock_pointer_fleets(0,f)-1);
	  int this_sel_model = selblock_models(selblock_pointer_fleets(0,f)-1);
	  temp_selAA = get_selAA_from_selAL(this_selAL, 0, this_sel_model, catch_phi_mat); // Original: should use catch_phi_mat
    for(int a = 0; a < n_ages; a++)
    {
      FAA(0,f,a) = F(0,f) * temp_selAA(a);
      log_FAA(0,f,a) = log(FAA(0,f,a));
      FAA_tot(0,a) = FAA_tot(0,a) + FAA(0,f,a);
    }
    for(int y = 1; y < n_years_model; y++)
    {
      log_F(y,f) = log_F(y-1,f) + F_devs(y-1,f);
      F(y,f) = exp(log_F(y,f));
	  // It transforms selAL to selAA if required using phi_matrix at fracyr.
	  matrix<Type> this_selAL = selAL(selblock_pointer_fleets(y,f)-1);
	  int this_sel_model = selblock_models(selblock_pointer_fleets(y,f)-1);
	  temp_selAA = get_selAA_from_selAL(this_selAL, y, this_sel_model, catch_phi_mat); // Original: should use catch_phi_mat
      for(int a = 0; a < n_ages; a++)
      {
        FAA(y,f,a) = F(y,f) * temp_selAA(a);
        log_FAA(y,f,a) = log(FAA(y,f,a));
        FAA_tot(y,a) = FAA_tot(y,a) + FAA(y,f,a);
      }
    }
  }
  // REPORT(temp_selAA);
  // Total mortality, Z = F + M (non-projection years only)
  for(int y = 0; y < n_years_model; y++) ZAA.row(y) = FAA_tot.row(y) + MAA.row(y);

  // ---------------------------------------------------------------------------------
  // Set up population model
  // Year 1 initialize population
  SSB.setZero();
  matrix<Type> NAA(n_years_model + n_years_proj,n_ages);
  NAA.setZero();
  matrix<Type> pred_NAA(n_years_model + n_years_proj,n_ages);
  pred_NAA.setZero();
  ssb_phi_mat = phi_matrix(waa_pointer_ssb-1);
  
  for(int a = 0; a < n_ages; a++)
  {
    if(N1_model == 0) NAA(0,a) = exp(log_N1_pars(a));
    else
    {
      if(a==0) NAA(0,0) = exp(log_N1_pars(0));
      else
      {
        if(a == n_ages-1) NAA(0,a) = NAA(0,a-1)/(1.0 + exp(-MAA(0,a) - exp(log_N1_pars(1)) * FAA_tot(0,a)/FAA_tot(0,which_F_age(0)-1)));
        else NAA(0,a) = NAA(0,a-1)* exp(-MAA(0,a) -  exp(log_N1_pars(1)) * FAA_tot(0,a)/FAA_tot(0,which_F_age(0)-1));
      }
    }
    // Calculate SSB using maturity at age:
	SSB(0) += NAA(0,a) * pred_waa(waa_pointer_ssb-1,0,a) * mat_at_age(0,a) * exp(-ZAA(0,a)*fracyr_SSB(0)); 
    pred_NAA(0,a) = NAA(0,a);
  }

  // get SPR0
  vector<Type> M(n_ages), sel(n_ages), mat(n_ages), waassb(n_ages), log_SPR0(n_years_model + n_years_proj);
  matrix<Type> waacatch(n_fleets,n_ages);
  int na = n_years_model + n_years_proj;
  vector<Type> log_SR_a(na), log_SR_b(na), SR_h(na), SR_R0(na);
  for(int y = 0; y < n_years_model + n_years_proj; y++)
  {
    for(int a = 0; a < n_ages; a++)
    {
      M(a) = MAA(y,a);
      waassb(a) = pred_waa(waa_pointer_ssb-1,y,a);
      mat(a) = mat_at_age(y,a);
    }
    log_SPR0(y) = log(get_SPR_0(M, mat, waassb, fracyr_SSB(y)));
  }
  REPORT(log_SPR0);

  // calculate stock-recruit parameters (steepness, R0, a, b)
  if(recruit_model > 2) //BH or Ricker SR
  {
    vector<Type> SR_h_tf(SR_h.size()); //different transformations for BH and Ricker
    if(recruit_model == 3) //BH stock recruit
    {
      if(use_steepness == 1)
      {
        SR_h.fill(0.2 + 0.8/(1+exp(-mean_rec_pars(0)))); //SR_a * SPR0/(4.0 + SR_a*SPR0);
        SR_R0.fill(exp(mean_rec_pars(1))); //(SR_a - 1/SPR0) / SR_b;
        log_SR_a = log(4 * SR_h/(exp(log_SPR0)*(1 - SR_h)));
        log_SR_b = log((5*SR_h - 1)/((1-SR_h)*SR_R0*exp(log_SPR0)));
      }
      else
      {
        log_SR_a.fill(mean_rec_pars(0));
        log_SR_b.fill(mean_rec_pars(1));
      }
      for(int i=0; i < n_Ecov; i++){
        if(Ecov_where(i,0) == 1){ // if ecov i affects recruitment
          for(int y = 0; y < n_years_model + n_years_proj; y++)
          {
            // (1) "controlling" = dens-indep mortality or (4) "masking" = metabolic/growth (decreases dR/dS)
            if((Ecov_how(i) == 1) | (Ecov_how(i) == 4))
            {
              log_SR_a(y) += Ecov_lm(i,0,y,0);
            }
            // (2) "limiting" = carrying capacity or (4) "masking" = metabolic/growth (decreases dR/dS)
            if((Ecov_how(i) == 2) | (Ecov_how(i) == 4))
            {
              log_SR_b(y) += Ecov_lm(i,0,y,0);
            }
          }
        }
      }
      if(use_steepness != 1)
      {
        SR_h = exp(log_SR_a) * exp(log_SPR0)/(4.0 + exp(log_SR_a + log_SPR0));
        SR_R0 = (exp(log_SR_a) - 1/exp(log_SPR0)) / exp(log_SR_b);
      }
      SR_h_tf = log(SR_h - 0.2) - log(1 - SR_h);
    }
    if(recruit_model>3) //Ricker stock recruit
    {
      if(use_steepness == 1)
      {
        SR_h.fill(0.2 + exp(mean_rec_pars(0)));
        SR_R0.fill(exp(mean_rec_pars(1)));
        log_SR_a = 1.25*log(5*SR_h) - log_SPR0;
        log_SR_b = log(1.25*log(5*SR_h)/(SR_R0*exp(log_SPR0)));
      }
      else
      {
        log_SR_a.fill(mean_rec_pars(0));
        log_SR_b.fill(mean_rec_pars(1));
      }
      for(int i=0; i < n_Ecov; i++){
        if(Ecov_where(i,0) == 1){ // if ecov i affects recruitment
          for(int y = 0; y < n_years_model + n_years_proj; y++)
          {
            if(Ecov_how(i) == 1) // "controlling" = dens-indep mortality
            {
              log_SR_a(y) += Ecov_lm(i,0,y,0);
            }
            if(Ecov_how(i) == 4) // "masking" = metabolic/growth (decreases dR/dS)
            { //NB: this is not identical to Iles and Beverton (1998), but their definition can give negative values of "b"
              log_SR_b(y) += 1.0 + Ecov_lm(i,0,y,0);
            }
          }
        }
      }
      if(use_steepness != 1)
      {
        SR_h = 0.2 * exp(0.8*log(exp(log_SR_a) * exp(log_SPR0)));
        SR_R0 = log(exp(log_SR_a + log_SPR0))/(exp(log_SR_b + log_SPR0));
      }
      SR_h_tf = log(SR_h - 0.2);
    }
    vector<Type> log_SR_R0 = log(SR_R0);
    if(do_post_samp.sum()==0){
      ADREPORT(log_SR_a);
      ADREPORT(log_SR_b);
      ADREPORT(SR_h_tf);
      ADREPORT(log_SR_R0);
    }
    REPORT(log_SR_a);
    REPORT(log_SR_b);
    REPORT(SR_h_tf);
    REPORT(log_SR_R0);
  }
  
  // ---------------------------------------------------------------------------------
  // Population model (get NAA, numbers-at-age, LAA, for all years)
  array<Type> NAA_devs(n_years_model+n_years_proj-1, n_ages);
  NAA_devs.setZero();

  for(int y = 1; y < n_years_model + n_years_proj; y++)
  {
    pred_NAA.row(y) = get_pred_NAA_y(y, recruit_model, mean_rec_pars, SSB, NAA, log_SR_a, 
      log_SR_b, Ecov_where, Ecov_how, Ecov_lm, ZAA);
    
    // calculate NAA
    if(n_NAA_sigma > 1){
      // all NAA are estimated (random effects)
      for(int a = 0; a < n_ages; a++) NAA(y,a) = exp(log_NAA(y-1,a));
      // calculate mean-0 deviations of log NAA (possibly bias-corrected)
      for(int a = 0; a < n_ages; a++) NAA_devs(y-1,a) = log_NAA(y-1,a) - log(pred_NAA(y,a));
    } else { // only recruitment estimated (either fixed or random effects)
      for(int a = 1; a < n_ages; a++) NAA(y,a) = pred_NAA(y,a); // for ages > 1 survival is deterministic 
      if((n_NAA_sigma == 0) && (y > n_years_model-1)){  //recruit FE, but recruit RE in projection years
        NAA(y,0) = exp(logR_proj(y-n_years_model)); // SCAA recruit in projections use diff object (random effect)
        for(int a = 1; a < n_ages; a++) NAA_devs(y-1,a) = log(NAA(y,a)) - log(pred_NAA(y,a));
        NAA_devs(y-1,0) = logR_proj(y-n_years_model) - log(pred_NAA(y,0));
      } else { //recruit RE, or (recruit FE and not projection year)
        NAA(y,0) = exp(log_NAA(y-1,0));
        for(int a = 0; a < n_ages; a++) NAA_devs(y-1,a) = log(NAA(y,a)) - log(pred_NAA(y,a));
      }
    }
        
    // calculate F and Z in projection years, here bc need NAA(y) if using F from catch
    if(do_proj == 1){ // now need FAA by fleet for projections, use total of average FAA by fleet over avg.yrs
      // get selectivity using average over avg.yrs
      if(y > n_years_model-1){
		waacatch = get_waacatch_y(pred_waa, y, n_ages, waa_pointer_fleets);
		waassb = get_waa_y(pred_waa, y, n_ages, waa_pointer_ssb);
        //n_fleets x n_ages: projected full F is sum of (means across years at age) across fleets 
        matrix<Type> FAA_proj = get_F_proj(y, n_fleets, proj_F_opt, FAA, NAA, MAA, mat_at_age, waacatch, waassb, fracyr_SSB, 
          log_SPR0, avg_years_ind, n_years_model, which_F_age, percentSPR, proj_Fcatch, percentFXSPR, F_proj_init(y-n_years_model), 
          log_SR_a, log_SR_b, recruit_model, percentFMSY);
        FAA_tot.row(y) = FAA_proj.colwise().sum();
        for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) FAA(y,f,a) = FAA_proj(f,a);
        //FAA_tot.row(y) = get_F_proj(y, proj_F_opt, FAA_tot, NAA, MAA, mat_at_age, waacatch, waassb, fracyr_SSB, log_SPR0, avg_years_ind, n_years_model,
        // which_F_age, percentSPR, proj_Fcatch);
        ZAA.row(y) = FAA_tot.row(y) + MAA.row(y);
      }
    } // end proj F
	// Calculate SSB:
	SSB(y) = get_SSB(NAA,ZAA,pred_waa,mat_at_age,y,waa_pointer_ssb,fracyr_SSB);
  } // end pop model loop

  // --------------------------------------------------------------------------
  // NAA random effects
  Type nll_NAA = Type(0);
  Type NAA_rho_a = rho_trans(trans_NAA_rho(0));
  Type NAA_rho_y = rho_trans(trans_NAA_rho(1));
  vector<Type> NAA_sigma = exp(log_NAA_sigma);
  vector<Type> sigma_a_sig(n_ages);
  sigma_a_sig.setZero();
  if(n_NAA_sigma == 1){
    sigma_a_sig(0) = NAA_sigma(0) / pow((1-pow(NAA_rho_y,2)),0.5);
  }
  if(n_NAA_sigma > 1){
    for(int a=0; a<n_ages; a++) sigma_a_sig(a) = NAA_sigma(NAA_sigma_pointers(a)-1) / pow((1-pow(NAA_rho_y,2))*(1-pow(NAA_rho_a,2)),0.5);
  }
  if(n_NAA_sigma > 0){
    if(do_post_samp.sum()==0){
      ADREPORT(NAA_sigma);
      ADREPORT(NAA_rho_a);
      ADREPORT(NAA_rho_y);
    }
    if(do_post_samp(0) == 1){
      ADREPORT(log_NAA);
    }
  }

  // likelihood of NAA deviations
  if((n_NAA_sigma == 0) && (do_proj == 1)){ // SCAA treats recruitment in proj years as random effects with fixed mean, SD
    for(int y = 0; y < n_years_proj; y++){
      nll_NAA -= dnorm(logR_proj(y), logR_mean, logR_sd, 1);
    }
    SIMULATE if((simulate_state(0) == 1) & (simulate_period(1) == 1)){ // proj years only
      for(int y = 0; y < n_years_proj; y++){
        logR_proj(y) = rnorm(logR_mean, logR_sd);
        NAA_devs(y+n_years_model-1,0) = logR_proj(y) - log(pred_NAA(y+n_years_model,0));
      }
    }
    REPORT(logR_proj);
  }
  if(n_NAA_sigma == 1){
    if(bias_correct_pe == 1) NAA_devs.col(0) += 0.5*pow(sigma_a_sig(0),2); //make sure this is ok when just recruitment is random.
    nll_NAA += SCALE(AR1(NAA_rho_y),sigma_a_sig(0))(NAA_devs.col(0));
    SIMULATE if(simulate_state(0) == 1) {
      vector<Type> NAAdevs0 = NAA_devs.col(0);
      AR1(NAA_rho_y).simulate(NAAdevs0); // sigma = 1, scale below
      NAAdevs0 = sigma_a_sig(0) * NAAdevs0;
      if(bias_correct_pe == 1) NAAdevs0 -= 0.5*pow(sigma_a_sig(0),2);
      for(int y = 0; y < n_years_model + n_years_proj - 1; y++){
        if(((simulate_period(0) == 1) & (y < n_years_model - 1)) | ((simulate_period(1) == 1) & (y > n_years_model - 2))){
          NAA_devs(y,0) = NAAdevs0(y);
        }
      }
    }
  }
  if(n_NAA_sigma > 1){
    if(bias_correct_pe == 1) for(int a = 0; a < n_ages; a++) NAA_devs.col(a) += 0.5*pow(sigma_a_sig(a),2);
    nll_NAA += SEPARABLE(VECSCALE(AR1(NAA_rho_a), sigma_a_sig),AR1(NAA_rho_y))(NAA_devs);
    SIMULATE if(simulate_state(0) == 1) {
      array<Type> NAAdevs = NAA_devs;
      SEPARABLE(VECSCALE(AR1(NAA_rho_a), sigma_a_sig),AR1(NAA_rho_y)).simulate(NAAdevs); // scaled here
      if(bias_correct_pe == 1) for(int a = 0; a < n_ages; a++) NAAdevs.col(a) -= 0.5*pow(sigma_a_sig(a),2);
      for(int y = 0; y < n_years_model + n_years_proj - 1; y++){
        if(((simulate_period(0) == 1) & (y < n_years_model - 1)) | ((simulate_period(1) == 1) & (y > n_years_model - 2))){
          for(int a = 0; a < n_ages; a++) NAA_devs(y,a) = NAAdevs(y,a);
        }
      }
    }
  }
  if((n_NAA_sigma > 0) | (do_proj == 1)) SIMULATE if(simulate_state(0) == 1){ // if n_NAA_sigma = 0 (SCAA), recruitment now random effects in projections
  	out_phi_mat = phi_matrix(waa_pointer_ssb-1);
	matrix<Type> sims = sim_pop(NAA_devs, recruit_model, mean_rec_pars, SSB,
	  NAA, log_SR_a, log_SR_b, Ecov_where, Ecov_how, Ecov_lm, 
      n_NAA_sigma, do_proj, proj_F_opt, FAA, FAA_tot, MAA, mat_at_age, pred_waa, waa_pointer_totcatch, waa_pointer_ssb, fracyr_SSB, log_SPR0, 
      avg_years_ind, n_years_model, n_fleets, which_F_age, percentSPR, proj_Fcatch, percentFXSPR, F_proj_init, percentFMSY);
    SSB = sims.col(sims.cols()-1);
    for(int a = 0; a < n_ages; a++) 
    {
      NAA.col(a) = sims.col(a);
      pred_NAA.col(a) = sims.col(a+n_ages);
      for(int y = 1;y < n_years_model + n_years_proj; y++) {
        log_NAA(y-1,a) = log(NAA(y,a));
        for(int f = 0; f < n_fleets; f++){
          FAA(y,f,a) = sims(y,2*n_ages + f*n_ages + a);
        }
      }
    }
    for(int y = 1;y < n_years_model + n_years_proj; y++) for(int a = 0; a < n_ages; a++){
      Type tot = 0;
      for(int f = 0; f < n_fleets; f++) tot += FAA(y,f,a);
      FAA_tot(y,a) = tot;
    }
    ZAA = MAA + FAA_tot;
    REPORT(sims);
    REPORT(log_NAA);
    REPORT(log_NAA_sigma);
    REPORT(trans_NAA_rho);
  }
  REPORT(NAA_devs);
  REPORT(nll_NAA);
  nll += nll_NAA;

  // ----------------------------------------------------------------
  // Catch data likelihood
  matrix<Type> nll_agg_catch(n_years_model,n_fleets), nll_catch_acomp(n_years_model,n_fleets), nll_catch_lcomp(n_years_model,n_fleets), agg_catch_proj(n_years_proj,n_fleets);
  array<Type> nll_catch_caal(n_years_model,n_fleets,n_lengths);
  array<Type> catch_paa_proj(n_fleets, n_years_proj, n_ages);
  array<Type> catch_pal_proj(n_fleets, n_years_proj, n_lengths);
  array<Type> catch_caal_proj(n_fleets,n_years_proj,n_lengths, n_ages); 
  matrix<Type> tmp_aging(n_ages,n_ages);
  vector<Type> tmp_agecomps(n_ages);
  nll_agg_catch.setZero();
  nll_catch_acomp.setZero();
  nll_catch_lcomp.setZero(); 
  nll_catch_caal.setZero();
  vector<Type> lsum(n_lengths);
  vector<Type> asum(n_ages);
  for(int y = 0; y < n_years_model+n_years_proj; y++)
  {
    //for now just use uncertainty from last year of catch
    int usey = y;
    if(y > n_years_model-1) usey = n_years_model-1;
    //int acomp_par_count = 0;
	for(int f = 0; f < n_fleets; f++)
    {
	  lsum.setZero();
	  asum.setZero();
      pred_catch(y,f) = 0.0;
      for(int a = 0; a < n_ages; a++){
		for(int l = 0; l < n_lengths; l++) { 
			// for model years, we can use either age or len selex: F vector goes until n_model_years
			if(y < n_years_model) pred_CAAL(y,f,l,a) = selAA(selblock_pointer_fleets(usey,f)-1)(usey,a)*selLL(selblock_pointer_fleets(usey,f)-1)(usey,l)*catch_phi_mat(l,a,y)*NAA(y,a)*F(y,f)*(1-exp(-ZAA(y,a)))/ZAA(y,a);
			// for projection years, we can only use selectivity at age as calculated above. TODO: implement len-selex for projection years.
			if(y > n_years_model-1) pred_CAAL(y,f,l,a) = catch_phi_mat(l,a,y)*NAA(y,a)*FAA(y,f,a)*(1-exp(-ZAA(y,a)))/ZAA(y,a);
			lsum(l) += pred_CAAL(y,f,l,a);
			asum(a) += pred_CAAL(y,f,l,a);
		}
	  }		
	  for(int l = 0; l < n_lengths; l++) pred_CAL(y,f,l) = lsum(l); // predicted catch-at-length
	  for(int a = 0; a < n_ages; a++) pred_CAA(y,f,a) = asum(a); // predicted catch-at-age (numbers)

	  // Calculate agg catch (biomass)
	  // When using L-W parameters: use w-at-len
	  //if(weight_model == 2) {
		//  for(int l = 0; l < n_lengths; l++) pred_catch(y,f) += watl(y,l)*pred_CAL(y,f,l); 
	  //}
	  // When using emp WAA or non parametric WAA: use w-at-age
	  //if((weight_model == 1) | (weight_model == 3)) {
	  for(int a = 0; a < n_ages; a++) pred_catch(y,f) += pred_waa(waa_pointer_fleets(f)-1,y,a) * pred_CAA(y,f,a); // biomass
	  //}

	  pred_log_catch(y,f) = log(pred_catch(y,f));
      Type sig = agg_catch_sigma(usey,f)*exp(log_catch_sig_scale(f));
      if(bias_correct_oe == 1) pred_log_catch(y,f) -= 0.5*exp(2*log(sig));
      if(y < n_years_model) if(use_agg_catch(y,f) == 1){
        nll_agg_catch(y,f) -= keep(keep_C(y,f)) * dnorm(obsvec(keep_C(y,f)), pred_log_catch(y,f), sig,1);
        nll_agg_catch(y,f) -= keep.cdf_lower(keep_C(y,f)) * log(squeeze(pnorm(obsvec(keep_C(y,f)), pred_log_catch(y,f), sig)));
        nll_agg_catch(y,f) -= keep.cdf_upper(keep_C(y,f)) * log(1.0 - squeeze(pnorm(obsvec(keep_C(y,f)), pred_log_catch(y,f), sig)));
      }
      SIMULATE if(simulate_data(0) == 1){
        if((simulate_period(0) == 1) & (y < n_years_model)) {
          agg_catch(y,f) = exp(rnorm(pred_log_catch(y,f), sig));
          if(use_agg_catch(y,f) == 1) obsvec(keep_C(y,f)) = log(agg_catch(y,f));
        }
        if((simulate_period(1) == 1) & (y > n_years_model - 1)) {
          agg_catch_proj(y-n_years_model,f) = exp(rnorm(pred_log_catch(y,f), sig));
        }
      }
	  
	  // catch PAA nll
      if(any_fleet_age_comp(f) == 1){					
        vector<Type> paa_obs_y(n_ages);
        paa_obs_y.setZero();
		if(use_catch_aging_error(f) == 1) { // use aging error
			for(int a = 0; a < n_ages; a++){
				for(int a2 = 0; a2 < n_ages; a2++) tmp_aging(a2,a) = pred_CAA(y,f,a)*catch_aging_error(f,a2,a);
			}
			tmp_agecomps = tmp_aging.rowwise().sum(); 
			for(int a = 0; a < n_ages; a++) {
				pred_catch_paa(y,f,a) = tmp_agecomps(a)/asum.sum(); // this object will contain the paa with aging error
				t_pred_paa(a) = pred_catch_paa(y,f,a);
			}
		} else { // not use aging error
			for(int a = 0; a < n_ages; a++){
				pred_catch_paa(y,f,a) = pred_CAA(y,f,a)/asum.sum();
				t_pred_paa(a) = pred_catch_paa(y,f,a);
			}
		}
			
        if(y < n_years_model) if(use_catch_paa(y,f) == 1) {
          for(int a = 0; a < n_ages; a++) paa_obs_y(a) = catch_paa(f,y,a);
          //NB: indexing in obsvec MUST be: keep_Cpaa(i,y,0),...,keep_Cpaa(i,y,0) + keep_Cpaa(i,y,1) - 1
          //keep_Cpaa(i,y,0) is first val, keep_Cpaa(i,y,1) is the length of the vector
          vector<Type> tf_paa_obs = obsvec.segment(keep_Cpaa(f,y,0), keep_Cpaa(f,y,1));
          vector<int> ages_obs_y = agesvec.segment(keep_Cpaa(f,y,0), keep_Cpaa(f,y,1));
          nll_catch_acomp(y,f) -= get_acomp_ll(tf_paa_obs, t_pred_paa, catch_Neff(y,f), ages_obs_y, age_comp_model_fleets(f), 
            vector<Type>(catch_paa_pars.row(f)), keep.segment(keep_Cpaa(f,y,0),keep_Cpaa(f,y,1)), do_osa, paa_obs_y);
        }
        SIMULATE if(simulate_data(0) == 1) if(use_catch_paa(usey,f) == 1){
          if((simulate_period(0) == 1) & (y < n_years_model)) //model years
          {
            for(int a = 0; a < n_ages; a++) paa_obs_y(a) = catch_paa(f,y,a);
            vector<int> ages_obs_y = agesvec.segment(keep_Cpaa(f,y,0), keep_Cpaa(f,y,1));
            vector<Type> tf_paa_obs = sim_acomp(t_pred_paa, catch_Neff(y,f), ages_obs_y, age_comp_model_fleets(f), vector<Type>(catch_paa_pars.row(f)));
            obsvec.segment(keep_Cpaa(f,y,0),keep_Cpaa(f,y,1)) = tf_paa_obs;
            paa_obs_y = make_paa(tf_paa_obs, age_comp_model_fleets(f), ages_obs_y, paa_obs_y);
            for(int a = 0; a < n_ages; a++) catch_paa(f,y,a) = paa_obs_y(a);
          }
          if((simulate_period(1) == 1) & (y > n_years_model - 1)) //projection years
          {
            for(int a = 0; a < n_ages; a++) paa_obs_y(a) = 1/n_ages; //only needed for LN obs to tell where the last non-zero age class is. No zeros in projections.
            vector<int> ages_obs_y(n_ages);
            for(int a = 0; a < n_ages; a++) ages_obs_y(a) = a;
            vector<Type> tf_paa_obs = sim_acomp(t_pred_paa, catch_Neff(usey,f), ages_obs_y, age_comp_model_fleets(f), vector<Type>(catch_paa_pars.row(f)));
            paa_obs_y = make_paa(tf_paa_obs, age_comp_model_fleets(f), ages_obs_y, paa_obs_y);
            for(int a = 0; a < n_ages; a++) catch_paa_proj(f,y-n_years_model,a) = paa_obs_y(a);
          }
        }
      }
	  
	  // catch PAL nll
      if(any_fleet_len_comp(f) == 1){
	    vector<Type> pal_obs_y(n_lengths);
        pal_obs_y.setZero();
        for(int l = 0; l < n_lengths; l++){
          pred_catch_pal(y,f,l) = pred_CAL(y,f,l)/lsum.sum();
          t_pred_pal(l) = pred_catch_pal(y,f,l);
        }
        if(y < n_years_model) if(use_catch_pal(y,f) == 1) {
          for(int l = 0; l < n_lengths; l++) pal_obs_y(l) = catch_pal(f,y,l);
          //NB: indexing in obsvec MUST be: keep_Cpaa(i,y,0),...,keep_Cpaa(i,y,0) + keep_Cpaa(i,y,1) - 1
          //keep_Cpaa(i,y,0) is first val, keep_Cpaa(i,y,1) is the length of the vector
          vector<Type> tf_pal_obs = obsvec.segment(keep_Cpal(f,y,0), keep_Cpal(f,y,1));
          nll_catch_lcomp(y,f) -= get_lcomp_ll(tf_pal_obs, t_pred_pal, catch_NeffL(y,f), len_comp_model_fleets(f), 
            vector<Type>(catch_pal_pars.row(f)), keep.segment(keep_Cpal(f,y,0),keep_Cpal(f,y,1)), do_osa, pal_obs_y);
        }
        SIMULATE if(simulate_data(0) == 1) if(use_catch_pal(usey,f) == 1){
          if((simulate_period(0) == 1) & (y < n_years_model)) //model years
          {
            for(int l = 0; l < n_lengths; l++) pal_obs_y(l) = catch_pal(f,y,l);
            vector<Type> tf_pal_obs = sim_lcomp(t_pred_pal, catch_NeffL(y,f), len_comp_model_fleets(f), vector<Type>(catch_pal_pars.row(f)));
            obsvec.segment(keep_Cpal(f,y,0),keep_Cpal(f,y,1)) = tf_pal_obs;
            pal_obs_y = make_pal(tf_pal_obs, len_comp_model_fleets(f));
            for(int l = 0; l < n_lengths; l++) catch_pal(f,y,l) = pal_obs_y(l);
          }
          if((simulate_period(1) == 1) & (y > n_years_model - 1)) //projection years
          {
            for(int l = 0; l < n_lengths; l++) pal_obs_y(l) = 1/n_lengths; //only needed for LN obs to tell where the last non-zero age class is. No zeros in projections.
            vector<Type> tf_pal_obs = sim_lcomp(t_pred_pal, catch_NeffL(usey,f), len_comp_model_fleets(f), vector<Type>(catch_pal_pars.row(f)));
            pal_obs_y = make_pal(tf_pal_obs, len_comp_model_fleets(f));
            for(int l = 0; l < n_lengths; l++) catch_pal_proj(f,y-n_years_model,l) = pal_obs_y(l);
          }
        }
      }
	  
	  // CAAL nll 
	  if(any_fleet_caal(f) == 1){
		  for(int l = 0; l < n_lengths; l++) {
			vector<Type> paa_obs_y(n_ages);
			paa_obs_y.setZero();
			if(use_catch_aging_error(f) == 1) { // use aging error
				for(int a = 0; a < n_ages; a++){
					for(int a2 = 0; a2 < n_ages; a2++) tmp_aging(a2,a) = pred_CAAL(y,f,l,a)*catch_aging_error(f,a2,a);
				}
				tmp_agecomps = tmp_aging.rowwise().sum(); 
				for(int a = 0; a < n_ages; a++) {
					pred_catch_caal(y,f,l,a) = tmp_agecomps(a)/lsum(l); // this object will contain the paa with aging error
					t_pred_paa(a) = pred_catch_caal(y,f,l,a);
				}
			} else { // not use aging error
				for(int a = 0; a < n_ages; a++){
				  pred_catch_caal(y,f,l,a) = pred_CAAL(y,f,l,a)/lsum(l);
				  t_pred_paa(a) = pred_catch_caal(y,f,l,a);
				}
			}
			if(y < n_years_model) if(use_catch_caal(y,f,l) == 1) {
			  for(int a = 0; a < n_ages; a++) paa_obs_y(a) = catch_caal(f,y,l,a);
			  //NB: indexing in obsvec MUST be: keep_Cpaa(i,y,0),...,keep_Cpaa(i,y,0) + keep_Cpaa(i,y,1) - 1
			  //keep_Cpaa(i,y,0) is first val, keep_Cpaa(i,y,1) is the length of the vector
			  vector<Type> tf_paa_obs = obsvec.segment(keep_Ccaal(f,y,l,0), keep_Ccaal(f,y,l,1));
			  vector<int> ages_obs_y = agesvec.segment(keep_Ccaal(f,y,l,0), keep_Ccaal(f,y,l,1));
			  nll_catch_caal(y,f,l) -= get_acomp_ll(tf_paa_obs, t_pred_paa, catch_caal_Neff(y,f,l), ages_obs_y, age_comp_model_fleets(f), 
				vector<Type>(catch_paa_pars.row(f)), keep.segment(keep_Ccaal(f,y,l,0),keep_Ccaal(f,y,l,1)), do_osa, paa_obs_y);
			}
			SIMULATE if(simulate_data(0) == 1) if(use_catch_caal(usey,f,l) == 1){
			  if((simulate_period(0) == 1) & (y < n_years_model)) //model years
			  {
				for(int a = 0; a < n_ages; a++) paa_obs_y(a) = catch_caal(f,y,l,a);
				vector<int> ages_obs_y = agesvec.segment(keep_Ccaal(f,y,l,0), keep_Ccaal(f,y,l,1));
				vector<Type> tf_paa_obs = sim_acomp(t_pred_paa, catch_caal_Neff(y,f,l), ages_obs_y, age_comp_model_fleets(f), vector<Type>(catch_paa_pars.row(f)));
				obsvec.segment(keep_Ccaal(f,y,l,0),keep_Ccaal(f,y,l,1)) = tf_paa_obs;
				paa_obs_y = make_paa(tf_paa_obs, age_comp_model_fleets(f), ages_obs_y, paa_obs_y);
				for(int a = 0; a < n_ages; a++) catch_caal(f,y,l,a) = paa_obs_y(a);
			  }
			  if((simulate_period(1) == 1) & (y > n_years_model - 1)) //projection years
			  {
				for(int a = 0; a < n_ages; a++) paa_obs_y(a) = 1/n_ages; //only needed for LN obs to tell where the last non-zero age class is. No zeros in projections.
				vector<int> ages_obs_y(n_ages);
				for(int a = 0; a < n_ages; a++) ages_obs_y(a) = a;
				vector<Type> tf_paa_obs = sim_acomp(t_pred_paa, catch_caal_Neff(usey,f,l), ages_obs_y, age_comp_model_fleets(f), vector<Type>(catch_paa_pars.row(f)));
				paa_obs_y = make_paa(tf_paa_obs, age_comp_model_fleets(f), ages_obs_y, paa_obs_y);
				for(int a = 0; a < n_ages; a++) catch_caal_proj(f,y-n_years_model,l,a) = paa_obs_y(a);
			  }
			}
		  } // length loop
      }
	 
    }
  }
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(0) == 1) REPORT(agg_catch);
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(0) == 1) REPORT(catch_paa);
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(0) == 1) REPORT(catch_pal);
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(0) == 1) REPORT(catch_caal); 
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(1) == 1) REPORT(agg_catch_proj);
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(1) == 1) REPORT(catch_paa_proj);
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(1) == 1) REPORT(catch_pal_proj); 
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(1) == 1) REPORT(catch_caal_proj);
  REPORT(nll_agg_catch);
  nll += nll_agg_catch.sum();
  //see(nll);
  REPORT(nll_catch_acomp);
  REPORT(nll_catch_lcomp);
  REPORT(nll_catch_caal);
  nll += nll_catch_acomp.sum();
  nll += nll_catch_lcomp.sum(); 
  nll += nll_catch_caal.sum(); 
  //see(nll);

  // -----------------------------------------------------------------------------
  // Index/survey data likelihood
  matrix<Type> nll_agg_indices(n_years_catch,n_indices), nll_index_acomp(n_years_catch,n_indices), nll_index_lcomp(n_years_catch,n_indices), agg_indices_proj(n_years_proj, n_indices);
  array<Type> nll_index_caal(n_years_catch,n_indices,n_lengths);
  array<Type> index_paa_proj(n_indices,n_years_proj,n_ages);
  array<Type> index_pal_proj(n_indices,n_years_proj,n_lengths); 
  array<Type> index_caal_proj(n_indices,n_years_proj,n_lengths, n_ages); 
  nll_agg_indices.setZero();
  nll_index_acomp.setZero();
  nll_index_lcomp.setZero();
  nll_index_caal.setZero();
  pred_indices.setZero();
  //Type temp_indices = 0.0;
  for(int y = 0; y < n_years_model + n_years_proj; y++)
  {
    int usey = y;
    if(y > n_years_model - 1) usey = n_years_model -1; //some things only go up to n_years_model-1
    //int acomp_par_count = 0;
    for(int i = 0; i < n_indices; i++)
    {
	  out_phi_mat = phi_matrix(waa_pointer_indices(i)-1); // just do it once to save time
	  lsum.setZero();
	  asum.setZero();
	  //temp_indices = 0.0;
      for(int a = 0; a < n_ages; a++) {
		for(int l = 0; l < n_lengths; l++) { 
			// Q effect here (SS does it on the total abundance, but I think it is same thing):
			pred_IAAL(y,i,l,a) = q(y,i)*selAA(selblock_pointer_indices(usey,i)-1)(usey,a)*selLL(selblock_pointer_indices(usey,i)-1)(usey,l)*out_phi_mat(l,a,y)*NAA(y,a)*exp(-ZAA(y,a) * fracyr_indices(usey,i));
			lsum(l) += pred_IAAL(y,i,l,a);
			asum(a) += pred_IAAL(y,i,l,a);
		}
	  }		
	  for(int l = 0; l < n_lengths; l++) pred_IAL(y,i,l) = lsum(l); // predicted index-at-length (numbers)
	  for(int a = 0; a < n_ages; a++) {
		  if(units_index_paa(i) == 1) pred_IAA(y,i,a) = pred_waa(waa_pointer_indices(i)-1,y,a) * asum(a); // predicted index-at-age (biomass)
		  else pred_IAA(y,i,a) = asum(a); // predicted index-at-age (numbers)
	  }
	
	  // Calculate agg index (numbers or biomass). Here add Q:
	  // When using L-W parameters:
	  //if(weight_model == 2) {
		//  for(int l = 0; l < n_lengths; l++){
		//	// here include Q effect:
		//	if(units_indices(i) == 1) pred_indices(y,i) += q(y,i)*watl(y,l)*pred_IAL(y,i,l); // biomass
		//	else pred_indices(y,i) += q(y,i)*pred_IAL(y,i,l); // numbers
		//  }
	  //}
	  // When using emp WAA or non parametric WAA
	  //if((weight_model == 1) | (weight_model == 3)) {
	  for(int a = 0; a < n_ages; a++) {
		  if(units_indices(i) == 1) pred_indices(y,i) += pred_waa(waa_pointer_indices(i)-1,y,a) * asum(a); // biomass, use asum to make sure we are using abundance
		  else pred_indices(y,i) += asum(a); // numbers, use asum to make sure we are using abundance
	  }
	  //}

	  // calculate Index NLL:
      pred_log_indices(y,i) = log(pred_indices(y,i));
      Type sig = agg_index_sigma(usey,i)*exp(log_index_sig_scale(i));
      if(bias_correct_oe == 1) pred_log_indices(y,i) -= 0.5*exp(2*log(sig));
      if(y < n_years_model) if(use_indices(y,i) == 1) {
        nll_agg_indices(y,i) -= keep(keep_I(y,i)) * dnorm(obsvec(keep_I(y,i)), pred_log_indices(y,i), sig, 1);
        nll_agg_indices(y,i) -= keep.cdf_lower(keep_I(y,i)) * log(squeeze(pnorm(obsvec(keep_I(y,i)), pred_log_indices(y,i), sig)));
        nll_agg_indices(y,i) -= keep.cdf_upper(keep_I(y,i)) * log(1.0 - squeeze(pnorm(obsvec(keep_I(y,i)), pred_log_indices(y,i), sig)));
      }
      SIMULATE if(simulate_data(1) == 1){
        if((simulate_period(0) == 1) & (y < n_years_model)) {
          agg_indices(y,i) = exp(rnorm(pred_log_indices(y,i), sig));
          if(use_indices(y,i) == 1) obsvec(keep_I(y,i)) = log(agg_indices(y,i));
        }
        if((simulate_period(1) == 1) & (y > n_years_model - 1)) agg_indices_proj(y-n_years_model,i) = exp(rnorm(pred_log_indices(y,i), sig));
      }
      
	  // acomp nll
      if(any_index_age_comp(i) == 1)
      {
        vector<Type> paa_obs_y(n_ages);
        paa_obs_y.setZero();
		if(use_index_aging_error(i) == 1) { // use aging error
			for(int a = 0; a < n_ages; a++){
				for(int a2 = 0; a2 < n_ages; a2++) tmp_aging(a2,a) = pred_IAA(y,i,a)*index_aging_error(i,a2,a);
			}
			tmp_agecomps = tmp_aging.rowwise().sum(); 
			for(int a = 0; a < n_ages; a++) {
				pred_index_paa(y,i,a) = tmp_agecomps(a)/asum.sum(); // this object will contain the paa with aging error
				t_pred_paa(a) = pred_index_paa(y,i,a);
			}
		} else { // not use aging error
			for(int a = 0; a < n_ages; a++) {
			  pred_index_paa(y,i,a) = pred_IAA(y,i,a)/asum.sum();
			  t_pred_paa(a) = pred_index_paa(y,i,a);
			}
		}
        if(y < n_years_model) if(use_index_paa(y,i) == 1) {
          for(int a = 0; a < n_ages; a++) paa_obs_y(a) = index_paa(i,y,a);
          //NB: indexing in obsvec MUST be: keep_Ipaa(i,y,0),...,keep_Ipaa(i,y,0) + keep_Ipaa(i,y,1) - 1
          //keep_Ipaa(i,y,0) is first val, keep_Ipaa(i,y,1) is the length of the vector
          vector<Type> tf_paa_obs = obsvec.segment(keep_Ipaa(i,y,0), keep_Ipaa(i,y,1));
          vector<int> ages_obs_y = agesvec.segment(keep_Ipaa(i,y,0), keep_Ipaa(i,y,1));
          nll_index_acomp(y,i) -= get_acomp_ll(tf_paa_obs, t_pred_paa, index_Neff(y,i), ages_obs_y, age_comp_model_indices(i), 
            vector<Type>(index_paa_pars.row(i)), keep.segment(keep_Ipaa(i,y,0),keep_Ipaa(i,y,1)), do_osa, paa_obs_y);
        }
        SIMULATE if(simulate_data(1) == 1) if(use_index_paa(usey,i) == 1){
          if((simulate_period(0) == 1) & (y < n_years_model)) //model years
          {
            for(int a = 0; a < n_ages; a++) paa_obs_y(a) = index_paa(i,y,a);
            vector<int> ages_obs_y = agesvec.segment(keep_Ipaa(i,y,0), keep_Ipaa(i,y,1));
            vector<Type> tf_paa_obs = sim_acomp(t_pred_paa, index_Neff(usey,i), ages_obs_y, age_comp_model_indices(i), vector<Type>(index_paa_pars.row(i)));//acomp_pars);
            obsvec.segment(keep_Ipaa(i,y,0),keep_Ipaa(i,y,1)) = tf_paa_obs;
            paa_obs_y = make_paa(tf_paa_obs, age_comp_model_indices(i), ages_obs_y, t_paa);
            for(int a = 0; a < n_ages; a++) index_paa(i,y,a) = paa_obs_y(a);
          }
          if((simulate_period(1) == 1) & (y > n_years_model - 1)) //projection years
          {
            for(int a = 0; a < n_ages; a++) paa_obs_y(a) = 1/n_ages; //only needed for LN obs to tell where the last non-zero age class is. No zeros in projections.
            vector<int> ages_obs_y(n_ages);
            for(int a = 0; a < n_ages; a++) ages_obs_y(a) = a;
            vector<Type> tf_paa_obs = sim_acomp(t_pred_paa, index_Neff(usey,i), ages_obs_y, age_comp_model_indices(i), vector<Type>(index_paa_pars.row(i)));//acomp_pars);
            paa_obs_y = make_paa(tf_paa_obs, age_comp_model_indices(i), ages_obs_y, paa_obs_y);
            for(int a = 0; a < n_ages; a++) index_paa_proj(i,y-n_years_model,a) = paa_obs_y(a);
          }
        }
      }

	  // lcomp nll
      if(any_index_len_comp(i) == 1)
      {
	  	vector<Type> pal_obs_y(n_lengths);
        pal_obs_y.setZero();
        for(int l = 0; l < n_lengths; l++)
        {
          pred_index_pal(y,i,l) = pred_IAL(y,i,l)/lsum.sum();
          t_pred_pal(l) = pred_index_pal(y,i,l);
        }
        if(y < n_years_model) if(use_index_pal(y,i) == 1) {
          for(int l = 0; l < n_lengths; l++) pal_obs_y(l) = index_pal(i,y,l);
          //NB: indexing in obsvec MUST be: keep_Cpaa(i,y,0),...,keep_Cpaa(i,y,0) + keep_Cpaa(i,y,1) - 1
          //keep_Cpaa(i,y,0) is first val, keep_Cpaa(i,y,1) is the length of the vector
          vector<Type> tf_pal_obs = obsvec.segment(keep_Ipal(i,y,0), keep_Ipal(i,y,1));
          nll_index_lcomp(y,i) -= get_lcomp_ll(tf_pal_obs, t_pred_pal, index_NeffL(y,i), len_comp_model_indices(i), 
					vector<Type>(index_pal_pars.row(i)), keep.segment(keep_Ipal(i,y,0),keep_Ipal(i,y,1)), do_osa, pal_obs_y);
        }
        SIMULATE if(simulate_data(1) == 1) if(use_index_pal(usey,i) == 1){
          if((simulate_period(0) == 1) & (y < n_years_model)) //model years
          {
            for(int l = 0; l < n_lengths; l++) pal_obs_y(l) = index_pal(i,y,l);
            vector<Type> tf_pal_obs = sim_lcomp(t_pred_pal, index_NeffL(usey,i), len_comp_model_indices(i), vector<Type>(index_pal_pars.row(i)));
            obsvec.segment(keep_Ipal(i,y,0),keep_Ipal(i,y,1)) = tf_pal_obs;
            pal_obs_y = make_pal(tf_pal_obs, len_comp_model_indices(i));
            for(int l = 0; l < n_lengths; l++) index_pal(i,y,l) = pal_obs_y(l);
          }
          if((simulate_period(1) == 1) & (y > n_years_model - 1)) //projection years
          {
            for(int l = 0; l < n_lengths; l++) pal_obs_y(l) = 1/n_lengths; //only needed for LN obs to tell where the last non-zero age class is. No zeros in projections.
            vector<Type> tf_pal_obs = sim_lcomp(t_pred_pal, index_NeffL(usey,i), len_comp_model_indices(i), vector<Type>(index_pal_pars.row(i)));
            pal_obs_y = make_pal(tf_pal_obs, len_comp_model_indices(i));
            for(int l = 0; l < n_lengths; l++) index_pal_proj(i,y-n_years_model,l) = pal_obs_y(l);
          }
        }
      }	  
	  
	  // CAAL nll
	  if(any_index_caal(i) == 1)
      {
		for(int l = 0; l < n_lengths; l++){
			vector<Type> paa_obs_y(n_ages);
			paa_obs_y.setZero();
			if(use_index_aging_error(i) == 1) { // use aging error
				for(int a = 0; a < n_ages; a++){
					for(int a2 = 0; a2 < n_ages; a2++) tmp_aging(a2,a) = pred_IAAL(y,i,l,a)*index_aging_error(i,a2,a);
				}
				tmp_agecomps = tmp_aging.rowwise().sum(); 
				for(int a = 0; a < n_ages; a++) {
					pred_index_caal(y,i,l,a) = tmp_agecomps(a)/lsum(l); // this object will contain the paa with aging error
					t_pred_paa(a) = pred_index_caal(y,i,l,a);
				}
			} else { // not use aging error
				for(int a = 0; a < n_ages; a++){
					pred_index_caal(y,i,l,a) = pred_IAAL(y,i,l,a)/lsum(l);
					t_pred_paa(a) = pred_index_caal(y,i,l,a); 
				}
			}
			if(y < n_years_model) if(use_index_caal(y,i,l) == 1) {
			  for(int a = 0; a < n_ages; a++) paa_obs_y(a) = index_caal(i,y,l,a);
			  //NB: indexing in obsvec MUST be: keep_Ipaa(i,y,0),...,keep_Ipaa(i,y,0) + keep_Ipaa(i,y,1) - 1
			  //keep_Ipaa(i,y,0) is first val, keep_Ipaa(i,y,1) is the length of the vector
			  vector<Type> tf_paa_obs = obsvec.segment(keep_Icaal(i,y,l,0), keep_Icaal(i,y,l,1));
			  vector<int> ages_obs_y = agesvec.segment(keep_Icaal(i,y,l,0), keep_Icaal(i,y,l,1));
			  nll_index_caal(y,i,l) -= get_acomp_ll(tf_paa_obs, t_pred_paa, index_caal_Neff(y,i,l), ages_obs_y, age_comp_model_indices(i), 
				vector<Type>(index_paa_pars.row(i)), keep.segment(keep_Icaal(i,y,l,0),keep_Icaal(i,y,l,1)), do_osa, paa_obs_y);
			}
			SIMULATE if(simulate_data(1) == 1) if(use_index_caal(usey,i,l) == 1){
			  if((simulate_period(0) == 1) & (y < n_years_model)) //model years
			  {
				for(int a = 0; a < n_ages; a++) paa_obs_y(a) = index_caal(i,y,l,a);
				vector<int> ages_obs_y = agesvec.segment(keep_Icaal(i,y,l,0), keep_Icaal(i,y,l,1));
				vector<Type> tf_paa_obs = sim_acomp(t_pred_paa, index_caal_Neff(usey,i,l), ages_obs_y, age_comp_model_indices(i), vector<Type>(index_paa_pars.row(i)));//acomp_pars);
				obsvec.segment(keep_Icaal(i,y,l,0),keep_Icaal(i,y,l,1)) = tf_paa_obs;
				paa_obs_y = make_paa(tf_paa_obs, age_comp_model_indices(i), ages_obs_y, t_paa);
				for(int a = 0; a < n_ages; a++) index_caal(i,y,l,a) = paa_obs_y(a);
			  }
			  if((simulate_period(1) == 1) & (y > n_years_model - 1)) //projection years
			  {
				for(int a = 0; a < n_ages; a++) paa_obs_y(a) = 1/n_ages; //only needed for LN obs to tell where the last non-zero age class is. No zeros in projections.
				vector<int> ages_obs_y(n_ages);
				for(int a = 0; a < n_ages; a++) ages_obs_y(a) = a;
				vector<Type> tf_paa_obs = sim_acomp(t_pred_paa, index_caal_Neff(usey,i,l), ages_obs_y, age_comp_model_indices(i), vector<Type>(index_paa_pars.row(i)));//acomp_pars);
				paa_obs_y = make_paa(tf_paa_obs, age_comp_model_indices(i), ages_obs_y, paa_obs_y);
				for(int a = 0; a < n_ages; a++) index_caal_proj(i,y-n_years_model,l,a) = paa_obs_y(a);
			  }
			}
		}
      }

	  
    }
  }
  SIMULATE if(simulate_data(1) == 1) if(simulate_period(0) == 1) REPORT(agg_indices);
  SIMULATE if(simulate_data(1) == 1) if(simulate_period(0) == 1) REPORT(index_paa);
  SIMULATE if(simulate_data(1) == 1) if(simulate_period(0) == 1) REPORT(index_pal);
  SIMULATE if(simulate_data(1) == 1) if(simulate_period(0) == 1) REPORT(index_caal);
  SIMULATE if(simulate_data(1) == 1) if(simulate_period(1) == 1) REPORT(agg_indices_proj);
  SIMULATE if(simulate_data(1) == 1) if(simulate_period(1) == 1) REPORT(index_paa_proj);
  SIMULATE if(simulate_data(1) == 1) if(simulate_period(1) == 1) REPORT(index_pal_proj);
  SIMULATE if(simulate_data(1) == 1) if(simulate_period(1) == 1) REPORT(index_caal_proj);
  REPORT(nll_agg_indices);
  nll += nll_agg_indices.sum();
  //see(nll);
  REPORT(nll_index_acomp);
  REPORT(nll_index_lcomp);
  REPORT(nll_index_caal);
  nll += nll_index_acomp.sum();
  nll += nll_index_lcomp.sum();
  nll += nll_index_caal.sum();
  //see(nll);
  
  SIMULATE if(sum(simulate_data) > 0) REPORT(obsvec);
  // -------------------------------------------------------------------
  // Calculate catch in projection years
  if(do_proj == 1){
    //vector<Type> catch_proj(n_years_proj), log_catch_proj(n_years_proj);
    matrix<Type> catch_proj(n_years_proj, n_fleets), log_catch_proj(n_years_proj, n_fleets);
    array<Type> CAA_proj(n_fleets, n_years_proj, n_ages);
    catch_proj.setZero();
    for(int i = 0; i < n_years_proj; i++){
      int yi = i + n_years_model;
      for(int a = 0; a < n_ages; a++){
        waacatch = get_waacatch_y(pred_waa, yi, n_ages, waa_pointer_fleets);
        for(int f = 0; f < n_fleets; f++) {
          CAA_proj(f,i,a) =  NAA(yi,a) * FAA(yi,f,a) * (1 - exp(-ZAA(yi,a)))/ZAA(yi,a);
          catch_proj(i,f) += waacatch(f,a) * CAA_proj(f,i,a);
        }
        for(int f = 0; f < n_fleets; f++) log_catch_proj(i,f) = log(catch_proj(i,f) + Type(1.0e-15));
      }
    }
    REPORT(catch_proj);
    if(do_post_samp.sum()==0) ADREPORT(log_catch_proj);
  }

  // ------------------------------------------------------------------------------

  //////////////////////////////////////////
  //Still need to add in yearly vectors of biological inputs, make sure to calculate SR_a,SR_b vector or otherwise.
  //////////////////////////////////////////
  //calculate BRPs
  //First SPR-based proxies
  //Type percentSPR = 40;
  vector<Type> predR(pred_NAA.rows());
  if(XSPR_R_opt == 1) predR = NAA.col(0);
  if(XSPR_R_opt == 3) predR = pred_NAA.col(0);
  if(XSPR_R_opt == 2){
    vector<Type> predR_toavg = XSPR_R_avg_yrs.unaryExpr(NAA.col(0));
    predR.fill(predR_toavg.mean());
  }
  if(XSPR_R_opt == 4){
    vector<Type> predR_toavg = XSPR_R_avg_yrs.unaryExpr(pred_NAA.col(0));
    predR.fill(predR_toavg.mean());
  }
  matrix<Type> SPR_res = get_SPR_res(MAA, FAA, which_F_age, pred_waa, waa_pointer_ssb, waa_pointer_fleets, mat_at_age, percentSPR, predR, fracyr_SSB, log_SPR0, FXSPR_init);
  vector<Type> log_FXSPR = SPR_res.col(0);
  vector<Type> log_SSB_FXSPR = SPR_res.col(1);
  vector<Type> log_Y_FXSPR = SPR_res.col(2);
  vector<Type> log_SPR_FXSPR = SPR_res.col(3);
  vector<Type> log_YPR_FXSPR = SPR_res.col(4);
  matrix<Type> log_FXSPR_iter = SPR_res.block(0,5,n_years_model + n_years_proj,10);

  REPORT(log_FXSPR_iter);
  REPORT(log_FXSPR);
  REPORT(log_SSB_FXSPR);
  REPORT(log_Y_FXSPR);
  REPORT(log_SPR_FXSPR);
  REPORT(log_YPR_FXSPR);
  if(do_post_samp.sum()==0){
    ADREPORT(log_FXSPR);
    ADREPORT(log_SSB_FXSPR);
    ADREPORT(log_Y_FXSPR);
  }

  //static/avg year results
  vector<Type> SPR_res_static = get_static_SPR_res(MAA, FAA, which_F_age_static, pred_waa, waa_pointer_ssb, waa_pointer_fleets, mat_at_age, percentSPR, NAA, 
    fracyr_SSB, static_FXSPR_init, avg_years_ind, avg_years_ind, avg_years_ind, avg_years_ind, avg_years_ind, XSPR_R_avg_yrs);
  Type log_FXSPR_static = SPR_res_static(0);
  Type log_SSB_FXSPR_static = SPR_res_static(1);
  Type log_Y_FXSPR_static = SPR_res_static(2);
  Type log_SPR_FXSPR_static = SPR_res_static(3);
  Type log_YPR_FXSPR_static = SPR_res_static(4);
  Type log_SPR0_static = SPR_res_static(5);
  vector<Type> log_FXSPR_iter_static = SPR_res_static.segment(6,10);

  REPORT(log_SPR0_static);
  REPORT(log_FXSPR_iter_static);
  REPORT(log_FXSPR_static);
  REPORT(log_SSB_FXSPR_static);
  REPORT(log_Y_FXSPR_static);
  REPORT(log_SPR_FXSPR_static);
  REPORT(log_YPR_FXSPR_static);
  if(do_post_samp.sum()==0){
    ADREPORT(log_FXSPR_static);
    ADREPORT(log_SSB_FXSPR_static);
    ADREPORT(log_Y_FXSPR_static);
  }

  //If stock-recruit models
  if(recruit_model > 2) //Beverton-Holt or Ricker selected
  {
    int n = 10;
    vector<Type> log_FMSY(n_years_model + n_years_proj), log_FMSY_i(1);
    matrix<Type> log_FMSY_iter(n_years_model + n_years_proj,n);
    vector<Type> log_YPR_MSY(n_years_model + n_years_proj), log_SPR_MSY(n_years_model + n_years_proj), log_R_MSY(n_years_model + n_years_proj);
    vector<Type> waacatch_MSY(n_ages);
    Type SR_a, SR_b;
    for(int y = 0; y < n_years_model + n_years_proj; y++)
    {
      log_FMSY_iter(y,0) = log(FMSY_init(y)); //starting value
      for(int a = 0; a < n_ages; a++)
      {
        M(a) = MAA(y,a);
        sel(a) = FAA_tot(y,a)/FAA_tot(y,which_F_age(y)-1); //have to look at FAA_tot to see where max F is.
        waassb(a) = pred_waa(waa_pointer_ssb-1,y,a);
        waacatch_MSY(a) = pred_waa(waa_pointer_totcatch-1, y, a);
        mat(a) = mat_at_age(y,a);
      }
      SR_a = exp(log_SR_a(y));
      SR_b = exp(log_SR_b(y));
      if(recruit_model == 3) //Beverton-Holt selected
      {
        sr_yield<Type> sryield(SR_a, SR_b, M, sel, mat, waassb, waacatch_MSY,fracyr_SSB(y),0);
        for (int i=0; i<n-1; i++)
        {
          log_FMSY_i(0) = log_FMSY_iter(y,i);
          vector<Type> grad_sr_yield = autodiff::gradient(sryield,log_FMSY_i);
          matrix<Type> hess_sr_yield = autodiff::hessian(sryield,log_FMSY_i);
          log_FMSY_iter(y,i+1) = log_FMSY_iter(y,i) - grad_sr_yield(0)/hess_sr_yield(0,0);
        }
      }
      else //Ricker selected
      {
        sr_yield<Type> sryield(SR_a, SR_b, M, sel, mat, waassb, waacatch_MSY,fracyr_SSB(y),1);
        for (int i=0; i<n-1; i++)
        {
          log_FMSY_i(0) = log_FMSY_iter(y,i);
          vector<Type> grad_sr_yield = autodiff::gradient(sryield,log_FMSY_i);
          matrix<Type> hess_sr_yield = autodiff::hessian(sryield,log_FMSY_i);
          log_FMSY_iter(y,i+1) = log_FMSY_iter(y,i) - grad_sr_yield(0)/hess_sr_yield(0,0);
        }
      }
      log_FMSY(y) = log_FMSY_iter(y,n-1);
      log_SPR_MSY(y) = log(get_SPR(log_FMSY(y), M, sel, mat, waassb, fracyr_SSB(y)));
      log_YPR_MSY(y) = log(get_YPR(log_FMSY(y), M, sel, waacatch_MSY));
      if(recruit_model == 3) log_R_MSY(y) = log((SR_a - 1/exp(log_SPR_MSY(y))) / SR_b); //bh
      else log_R_MSY(y) = log(log(SR_a) + log_SPR_MSY(y)) - log(SR_b) - log_SPR_MSY(y); //ricker
    }
    vector<Type> log_SSB_MSY = log_R_MSY + log_SPR_MSY;
    vector<Type> log_MSY = log_R_MSY + log_YPR_MSY;

    if(do_post_samp.sum()==0){
      ADREPORT(log_FMSY);
      ADREPORT(log_SSB_MSY);
      ADREPORT(log_R_MSY);
      ADREPORT(log_MSY);
      ADREPORT(log_SPR_MSY);
      ADREPORT(log_YPR_MSY);
    }
    REPORT(log_FMSY);
    REPORT(log_FMSY_iter);
    REPORT(log_SSB_MSY);
    REPORT(log_R_MSY);
    REPORT(log_MSY);
    REPORT(log_SPR_MSY);
    REPORT(log_YPR_MSY);
  }

  matrix<Type> log_FAA_tot = log(FAA_tot.array());
  matrix<Type> log_index_resid(n_years_model, n_indices);
  log_index_resid.setZero();
  for(int y = 0; y < n_years_model; y++){
    for(int i = 0; i < n_indices; i++){
      if(use_indices(y,i) == 1) log_index_resid(y,i) = log(agg_indices(y,i)) - pred_log_indices(y,i);
    }
  }
  matrix<Type> log_catch_resid = log(agg_catch.block(0,0,n_years_model,n_fleets).array()) - pred_log_catch.block(0,0,n_years_model,n_fleets).array();
  vector<Type> log_SSB =  log(SSB);
  vector<Type> Fbar(n_years_model + n_years_proj);
  Fbar.setZero();
  int n_Fbar_ages = Fbar_ages.size();
  for(int y = 0; y < n_years_model + n_years_proj; y++) for(int a = 0; a < n_Fbar_ages; a++) Fbar(y) += FAA_tot(y,Fbar_ages(a)-1)/n_Fbar_ages;

  vector<Type> log_Fbar = log(Fbar);
  matrix<Type> log_NAA_rep = log(NAA.array());

  REPORT(NAA);
  REPORT(jan1LAA); 
  REPORT(SDAA);
  REPORT(LAA_par); 
  REPORT(WAA_par); 
  REPORT(pred_NAA);
  REPORT(SSB);
  REPORT(selAL);
  REPORT(selAA);
  REPORT(jan1_phi_mat); 
  REPORT(ssb_phi_mat);
  REPORT(catch_phi_mat);
  REPORT(mat_at_age); 
  REPORT(selLL);
  REPORT(MAA);
  REPORT(F);
  REPORT(FAA);
  REPORT(FAA_tot);
  REPORT(ZAA);
  REPORT(Fbar);
  REPORT(pred_catch); // baranov eq, not bias-corrected
  REPORT(pred_log_catch); // bias-corrected
  REPORT(pred_catch_paa);
  REPORT(pred_catch_pal);
  REPORT(pred_catch_caal);
  REPORT(pred_CAA);
  REPORT(pred_CAL);
  REPORT(pred_CAAL);
  REPORT(pred_indices);
  REPORT(pred_log_indices); // bias-corrected
  REPORT(pred_index_paa);
  REPORT(pred_index_pal);
  REPORT(pred_index_caal);
  REPORT(pred_IAA);
  REPORT(pred_IAL);
  REPORT(pred_IAAL);
  REPORT(Ecov_x);
  REPORT(Ecov_out);
  REPORT(Ecov_process_pars);
  REPORT(Ecov_re);
  REPORT(Ecov_beta);
  REPORT(mean_rec_pars);
  REPORT(Ecov_obs_sigma_par);
  REPORT(Ecov_obs_sigma);

  if(do_post_samp.sum()==0){
    ADREPORT(log_F);
    ADREPORT(log_FAA);
    ADREPORT(log_FAA_tot);
    ADREPORT(log_Fbar);
    ADREPORT(log_NAA_rep);
    ADREPORT(log_SSB);
    ADREPORT(log_index_resid);
    ADREPORT(log_catch_resid);
  }

  REPORT(nll);
  REPORT(nll_Ecov);
  REPORT(nll_Ecov_obs);

  if(Ecov_model.sum() > 0){
    matrix<Type> Ecov_resid = Ecov_obs.array() - Ecov_x.block(0,0,n_years_Ecov,n_Ecov).array();
    if(do_post_samp.sum()==0){
      ADREPORT(Ecov_x);
      ADREPORT(Ecov_resid);
    }
  }
  
  // DM linear (age comps):
  // for fisheries:
  vector<int> any_DM_fleets(n_fleets);
  any_DM_fleets.setZero();
  for(int f = 0; f < n_fleets; f++){
    if(age_comp_model_fleets(f) == 2) any_DM_fleets(f) = 1;
    if(age_comp_model_fleets(f) == 11) any_DM_fleets(f) = 2;
  }
  if(sum(any_DM_fleets)>0){
    matrix<Type> Neff_est_fleets(n_years_model, n_fleets);
    Neff_est_fleets.setZero();
    for(int f = 0; f < n_fleets; f++){
      if(any_DM_fleets(f) == 1) {
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) - log(N_eff) for normal D-M option, so CI's could be created from that SE estimate.
        for(int y = 0; y < n_years_model; y++) Neff_est_fleets(y,f) = 1 + (catch_Neff(y,f) -1)/(1 + catch_Neff(y,f) * exp(-catch_paa_pars(f,0)));
      }
      if(any_DM_fleets(f) == 2) {
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) for linear D-M option, so CI's could be created from that SE estimate.
        for(int y = 0; y < n_years_model; y++) Neff_est_fleets(y,f) = 1 + (catch_Neff(y,f) -1)/(1 + exp(-catch_paa_pars(f,0)));
      }
    }
    REPORT(Neff_est_fleets);
  }
  // for indices:
  vector<int> any_DM_indices(n_indices);
  any_DM_indices.setZero();
  for(int i = 0; i < n_indices; i++){
    if(age_comp_model_indices(i) == 2) any_DM_indices(i) = 1;
    if(age_comp_model_indices(i) == 11) any_DM_indices(i) = 2;
  }
  if(sum(any_DM_indices)>0){
    matrix<Type> Neff_est_indices(n_years_model, n_indices);
    Neff_est_indices.setZero();
    for(int i = 0; i < n_indices; i++){
      if(any_DM_fleets(i) == 1) {
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) - log(N_eff) for normal D-M option, so CI's could be created from that SE estimate.
        for(int y = 0; y < n_years_model; y++) Neff_est_indices(y,i) = 1 + (index_Neff(y,i) -1)/(1 + index_Neff(y,i)*exp(-index_paa_pars(i,0)));
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) for either D-M option, so CI's could be created from that SE estimate.
      }
      if(any_DM_fleets(i) == 2) {
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) for linear D-M option, so CI's could be created from that SE estimate.
        for(int y = 0; y < n_years_model; y++) Neff_est_indices(y,i) = 1 + (index_Neff(y,i) -1)/(1 + exp(-index_paa_pars(i,0)));
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) for either D-M option, so CI's could be created from that SE estimate.
      }
    }
    REPORT(Neff_est_indices);
  }
  
  // DM linear (len comps):
  // for fisheries:
  any_DM_fleets.setZero();
  for(int f = 0; f < n_fleets; f++){
    if(len_comp_model_fleets(f) == 2) any_DM_fleets(f) = 1;
    if(len_comp_model_fleets(f) == 3) any_DM_fleets(f) = 2;
  }
  if(sum(any_DM_fleets)>0){
    matrix<Type> NeffL_est_fleets(n_years_model, n_fleets);
    NeffL_est_fleets.setZero();
    for(int f = 0; f < n_fleets; f++){
      if(any_DM_fleets(f) == 1) {
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) - log(N_eff) for normal D-M option, so CI's could be created from that SE estimate.
        for(int y = 0; y < n_years_model; y++) NeffL_est_fleets(y,f) = 1 + (catch_NeffL(y,f) -1)/(1 + catch_NeffL(y,f) * exp(-catch_pal_pars(f,0)));
      }
      if(any_DM_fleets(f) == 2) {
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) for linear D-M option, so CI's could be created from that SE estimate.
        for(int y = 0; y < n_years_model; y++) NeffL_est_fleets(y,f) = 1 + (catch_NeffL(y,f) -1)/(1 + exp(-catch_pal_pars(f,0)));
      }
    }
    REPORT(NeffL_est_fleets);
  }
  // for indices:
  any_DM_indices.setZero();
  for(int i = 0; i < n_indices; i++){
    if(len_comp_model_indices(i) == 2) any_DM_indices(i) = 1;
    if(len_comp_model_indices(i) == 3) any_DM_indices(i) = 2;
  }
  if(sum(any_DM_indices)>0){
    matrix<Type> NeffL_est_indices(n_years_model, n_indices);
    NeffL_est_indices.setZero();
    for(int i = 0; i < n_indices; i++){
      if(any_DM_fleets(i) == 1) {
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) - log(N_eff) for normal D-M option, so CI's could be created from that SE estimate.
        for(int y = 0; y < n_years_model; y++) NeffL_est_indices(y,i) = 1 + (index_NeffL(y,i) -1)/(1 + index_NeffL(y,i)*exp(-index_pal_pars(i,0)));
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) for either D-M option, so CI's could be created from that SE estimate.
      }
      if(any_DM_fleets(i) == 2) {
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) for linear D-M option, so CI's could be created from that SE estimate.
        for(int y = 0; y < n_years_model; y++) NeffL_est_indices(y,i) = 1 + (index_NeffL(y,i) -1)/(1 + exp(-index_pal_pars(i,0)));
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) for either D-M option, so CI's could be created from that SE estimate.
      }
    }
    REPORT(NeffL_est_indices);
  }
  
  // TODO: do it for CAAL

  SIMULATE {
    REPORT(logit_q);
    REPORT(log_F1);
    REPORT(F_devs);
    REPORT(log_N1_pars);
    REPORT(catch_paa_pars);
    REPORT(catch_pal_pars);
    REPORT(index_paa_pars);
    REPORT(index_pal_pars);
    REPORT(log_b);
    REPORT(log_catch_sig_scale);
    REPORT(log_index_sig_scale);
  }

  return nll;
}


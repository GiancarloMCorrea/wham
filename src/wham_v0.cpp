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
  DATA_INTEGER(n_lengths); // NEWG
  DATA_VECTOR(lengths); // NEWG
  DATA_INTEGER(n_fleets);
  DATA_INTEGER(n_indices);
  DATA_INTEGER(n_selblocks);
  DATA_IVECTOR(selblock_models); // for each block: 1 = age-specific, 2 = logistic, 3 = double-logistic, 4 = logistic (declining)
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
  DATA_IVECTOR(len_comp_model_fleets); // NEWG
  DATA_IVECTOR(len_comp_model_indices); // NEWG
  DATA_VECTOR(fracyr_SSB);
  DATA_MATRIX(mature);
  DATA_IVECTOR(waa_pointer_fleets);
  DATA_INTEGER(waa_pointer_totcatch);
  DATA_IVECTOR(waa_pointer_indices);
  DATA_INTEGER(waa_pointer_ssb);
  DATA_INTEGER(waa_pointer_jan1);
  DATA_INTEGER(waa_type);
  DATA_IMATRIX(use_catch_waa);
  DATA_IMATRIX(use_index_waa);
  DATA_ARRAY(waa);
  DATA_ARRAY(waa_cv);
  DATA_MATRIX(agg_catch);
  DATA_IMATRIX(use_agg_catch);
  DATA_MATRIX(agg_catch_sigma);
  DATA_ARRAY(catch_paa); //n_fleets x n_years x n_ages
  DATA_IMATRIX(use_catch_paa);
  DATA_MATRIX(catch_Neff);
  
  DATA_ARRAY(catch_pal); //n_fleets x n_years x n_lengths // NEWG
  DATA_IMATRIX(use_catch_pal); // NEWG
  DATA_MATRIX(catch_NeffL); // NEWG
  DATA_ARRAY(catch_caal); //n_fleets x n_years x n_lengths // NEWG
  DATA_IMATRIX(use_catch_caal); // NEWG
  DATA_ARRAY(catch_caal_Neff); // NEWG
  
  DATA_IVECTOR(units_indices);
  DATA_MATRIX(fracyr_indices);
  DATA_MATRIX(agg_indices);
  DATA_IMATRIX(use_indices);
  DATA_MATRIX(agg_index_sigma);
  DATA_IVECTOR(units_index_paa);
  DATA_ARRAY(index_paa); //n_indices x n_years x n_ages
  DATA_IMATRIX(use_index_paa);
  DATA_MATRIX(index_Neff);
  
  DATA_IVECTOR(units_index_pal); // NEWG
  DATA_ARRAY(index_pal); //n_indices x n_years x n_lengths // NEWG
  DATA_IMATRIX(use_index_pal); // NEWG
  DATA_MATRIX(index_NeffL);// NEWG
  DATA_ARRAY(index_caal); //n_fleets x n_years x n_lengths // NEWG
  DATA_IMATRIX(use_index_caal); // NEWG
  DATA_ARRAY(index_caal_Neff); // NEWG
  
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
  
  DATA_IVECTOR(LAA_est); //  NEWG
  DATA_INTEGER(LAA_re_model); // NEWG

  DATA_INTEGER(n_growth_par); // NEWG TODO: change this parameter name
  DATA_INTEGER(growth_model); // 1: "vB-classic", 2: "vB-K_age" NEWG
  DATA_IVECTOR(growth_re_model); // 1 = none, 2 = IID, 3 = ar1_y NEWG
  DATA_ARRAY(phi_matrix_input); // 1 = none, 2 = IID, 3 = ar1_y NEWG
  DATA_INTEGER(phi_matrix_info); // 1 = none, 2 = IID, 3 = ar1_y NEWG

  DATA_INTEGER(n_LW_par); // NEWG TODO: change this parameter name
  DATA_IVECTOR(LW_re_model); // 1 = none, 2 = IID, 3 = ar1_y NEWG
  
  DATA_IVECTOR(which_F_age); //which age of F to use for full total F for msy/ypr calculations and projections (n_years_model + n_years_proj)
  DATA_INTEGER(use_steepness); // which parameterization to use for BH/Ricker S-R, if needed.
  DATA_INTEGER(bias_correct_pe); //bias correct lognormal process error?
  DATA_INTEGER(bias_correct_oe); //bias correct lognormal observation error?
  DATA_IVECTOR(Fbar_ages);
  DATA_IVECTOR(simulate_state); //vector (0/1) if 1 then state parameters (NAA, MAA, sel, Ecov, q, growth, LW) in that order) will be simulated.
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
  DATA_IARRAY(keep_Cpal); // NEWG
  DATA_IARRAY(keep_Ipal); // NEWG
  DATA_IARRAY(keep_Ccaal); // NEWG
  DATA_IARRAY(keep_Icaal); // NEWG
  DATA_IVECTOR(do_post_samp); //length = 5, whether to ADREPORT posterior residuals for NAA, M, selectivity, Ecov, q. 

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
  DATA_IMATRIX(Ecov_where); // n_Ecov x 2+n_indices. 0/1 values with columns corresponding to recruit, mortality, indices in that order
  DATA_INTEGER(Ecov_where_subindex); // n_Ecov x 2+n_indices. 0/1 values with columns corresponding to recruit, mortality, indices in that order
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
  DATA_IVECTOR(proj_GW_opt); // 1 = continue M_re (check for time-varying M_re on R side), 2 = average M (over avg_years_ind)
  DATA_IVECTOR(proj_LW_opt); // 1 = continue M_re (check for time-varying M_re on R side), 2 = average M (over avg_years_ind)
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
  PARAMETER_MATRIX(catch_pal_pars); //n_fleets x 3 // NEWG
  //PARAMETER_VECTOR(index_paa_pars);
  PARAMETER_MATRIX(index_paa_pars); //n_indices x 3
  PARAMETER_MATRIX(index_pal_pars); //n_indices x 3 // NEWG
  PARAMETER_VECTOR(M_a); // mean M-at-age, fixed effects, length = n_ages if M_model = 2 (age-specific), length = 1 if M_model = 1 (constant) or 3 (weight-at-age M)
  PARAMETER_ARRAY(M_re); // random effects for year- and age-varying M deviations from mean M_a), dim = n_years x n_M_a
  PARAMETER_VECTOR(M_repars); // parameters controlling M_re, length = 3 (sigma_M, rho_M_a, rho_M_y)
  
  // NEWG section: growth parameters:
  PARAMETER_MATRIX(growth_a); //  
  PARAMETER_ARRAY(growth_re); //    
  PARAMETER_MATRIX(growth_repars); //    

  PARAMETER_VECTOR(LAA_a); //  
  PARAMETER_ARRAY(LAA_re); //  
  PARAMETER_VECTOR(LAA_repars); //   

  PARAMETER_MATRIX(LW_a); //  could be 3 parameters or input_LAA
  PARAMETER_ARRAY(LW_re); // nyears x nages x npars_growth     
  PARAMETER_MATRIX(LW_repars); //  npars_growth x npars_re     
    
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
  vector<int> any_index_len_comp(n_indices);  // NEWG
  vector<int> any_fleet_len_comp(n_fleets);  // NEWG
  vector<int> any_index_caal(n_indices);  // NEWG
  vector<int> any_fleet_caal(n_fleets);  // NEWG
  vector<Type> SSB(n_years_model + n_years_proj);
  matrix<Type> F(n_years_model,n_fleets);
  matrix<Type> log_F(n_years_model,n_fleets);
  
  array<Type> pred_CAA(n_years_model+n_years_proj,n_fleets,n_ages);
  array<Type> pred_catch_paa(n_years_model+n_years_proj,n_fleets,n_ages);
  array<Type> pred_CAL(n_years_model+n_years_proj,n_fleets,n_lengths); // NEWG
  array<Type> pred_CAAL(n_years_model+n_years_proj,n_fleets,n_lengths,n_ages); // NEWG
  array<Type> pred_catch_pal(n_years_model+n_years_proj,n_fleets,n_lengths); // NEWG
  array<Type> pred_catch_caal(n_years_model+n_years_proj,n_fleets,n_lengths,n_ages); // NEWG
  matrix<Type> pred_catch(n_years_model+n_years_proj,n_fleets);
  matrix<Type> pred_log_catch(n_years_model+n_years_proj,n_fleets);
  
  array<Type> pred_IAA(n_years_model+n_years_proj,n_indices,n_ages);
  array<Type> pred_IAAL(n_years_model+n_years_proj,n_indices,n_lengths,n_ages); // NEWG
  array<Type> pred_index_paa(n_years_model+n_years_proj,n_indices,n_ages);
  array<Type> pred_IAL(n_years_model+n_years_proj,n_indices,n_lengths); // NEWG
  array<Type> pred_index_pal(n_years_model+n_years_proj,n_indices,n_lengths); // NEWG
  array<Type> pred_index_caal(n_years_model+n_years_proj,n_indices,n_lengths,n_ages); // NEWG
  matrix<Type> pred_indices(n_years_model+n_years_proj,n_indices); // not bias corrected
  matrix<Type> pred_log_indices(n_years_model+n_years_proj,n_indices); // bias corrected  

  array<Type> pred_waa(waa.dim(0), n_years_model + n_years_proj, n_ages); //
  array<Type> phi_mat(n_years_model + n_years_proj,n_lengths,n_ages);// NEWG
  matrix<Type> LAA(n_years_model + n_years_proj,n_ages);// NEWG
  matrix<Type> fracyr_phi_mat(n_lengths,n_ages);// NEWG
  matrix<Type> SDAA(n_years_model + n_years_proj,n_ages);// NEWG
  array<Type> FAA(n_years_model+n_years_proj,n_fleets,n_ages);
  array<Type> log_FAA(n_years_model+n_years_proj,n_fleets,n_ages);
  matrix<Type> FAA_tot(n_years_model + n_years_proj,n_ages);
  matrix<Type> ZAA(n_years_model + n_years_proj,n_ages);
  array<Type> QAA(n_years_model+n_years_proj,n_indices,n_ages);
  vector<matrix<Type> > selAL(n_selblocks); // Could be selex at age or len
  vector<matrix<Type> > selAA(n_selblocks); // selAA(b)(y,a) gives selectivity by block, year, age; selAA(b) is matrix with dim = n_years x n_ages;
  matrix<Type> q(n_years_model+n_years_proj,n_indices);
  vector<Type> t_paa(n_ages); // used also for ALK
  vector<Type> t_pred_paa(n_ages); // used also for ALK
  vector<Type> t_pal(n_lengths); // NEWG
  vector<Type> t_pred_pal(n_lengths); // NEWG
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
  // NEWG section:
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
  // NEWG section ALK:
  for(int i = 0; i < n_indices; i++)
  {
    any_index_caal(i) = 0;
    for(int y = 0; y < n_years_indices; y++) if(use_index_caal(y,i) == 1) any_index_caal(i) = 1;
  }
  for(int i = 0; i < n_fleets; i++)
  {
    any_fleet_caal(i) = 0;
    for(int y = 0; y < n_years_catch; y++) if(use_catch_caal(y,i) == 1) any_fleet_caal(i) = 1;
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
  selAL = get_selectivity(n_years_model, n_ages, n_lengths, lengths, n_selblocks, selpars, selblock_models); // Get selectivity by block, age, year. This could be selex-at-len as well
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
  int n_effects = Ecov_beta.dim(0); // 2 + n_indices (recruitment, mortality and any catchabilities)
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

  // --------------------------------------------------------------------------
  // Calculate GROWTH --  NEWG section
  
  // Growth parameters re section
  Type nll_GW = Type(0);
  
  // For all growth parameters:
  vector<Type> sigma_GW(n_growth_par); // first RE parameter
  vector<Type> rho_GW_y(n_growth_par);  // second RE parameter

  for(int j = 0; j < n_growth_par; j++) {
  
  Type Sigma_GW = Type(0);

	  if((growth_re_model(j) == 2) | (growth_re_model(j) == 4)) { // Only for ii_y and Ar1_y. CHECK THIS LATER, ONLY FOR WORK AR_y so far I think
		
		sigma_GW(j) = exp(growth_repars(j,0));
        rho_GW_y(j) = rho_trans(growth_repars(j,1));  
		// likelihood of growth parameters deviations
			vector<Type> GWre0(n_years_model + n_years_proj);
			for(int y = 0; y < n_years_model + n_years_proj; y++) GWre0(y) = growth_re(y,0,j); 
			Sigma_GW = pow(pow(sigma_GW(j),2) / (1-pow(rho_GW_y(j),2)),0.5);
			nll_GW += SCALE(AR1(rho_GW_y(j)), Sigma_GW)(GWre0);
			SIMULATE if(simulate_state(5) == 1) if(sum(simulate_period) > 0) {
			  AR1(rho_GW_y(j)).simulate(GWre0);
			  for(int i = 0; i < GWre0.size(); i++) GWre0(i) = Sigma_GW * GWre0(i);
			  for(int y = 0; y < n_years_model + n_years_proj; y++){
				if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
				  for(int a = 0; a < n_ages; a++) growth_re(y,a,j) = GWre0(y);
				}
			  }
			}

		if(do_post_samp.sum()==0){
		  ADREPORT(sigma_GW);
		  ADREPORT(rho_GW_y);
		}
		
	  }
	  
	  if((growth_re_model(j) == 3) | (growth_re_model(j) == 5)) { // Only for ii_c and Ar1_c
		  
		sigma_GW(j) = exp(growth_repars(j,0));
        rho_GW_y(j) = rho_trans(growth_repars(j,1));  
		// likelihood of growth parameters deviations
			vector<Type> GWre0(n_years_model + n_years_proj + n_ages - 1);
			for(int i = 0; i < (n_ages - 1); i++) GWre0(i) = growth_re(0,n_ages - i - 1,j); // for cohorts at y = 0 except a = 0
			for(int i = (n_ages - 1); i < (n_years_model + n_years_proj + n_ages - 1); i++) GWre0(i) = growth_re(i-n_ages+1,0,j); // for cohorts y>=0 
			Sigma_GW = pow(pow(sigma_GW(j),2) / (1-pow(rho_GW_y(j),2)),0.5);
			nll_GW += SCALE(AR1(rho_GW_y(j)), Sigma_GW)(GWre0);
			SIMULATE if(simulate_state(5) == 1) if(sum(simulate_period) > 0) {
			  AR1(rho_GW_y(j)).simulate(GWre0);
			  for(int i = 0; i < GWre0.size(); i++) GWre0(i) = Sigma_GW * GWre0(i);
			  for(int y = 0; y < n_years_model + n_years_proj; y++){
				if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
				  for(int a = 0; a < n_ages; a++) {
					  if(y == 0) growth_re(y,n_ages-1-a,j) = GWre0(a); 
					  else growth_re(y,a,j) = GWre0(y-a+n_ages-1);
				  }
				}
			  }
			}

		if(do_post_samp.sum()==0){
		  ADREPORT(sigma_GW);
		  ADREPORT(rho_GW_y);
		}
		  
	  }
  
  }
  
  REPORT(nll_GW);
  nll += nll_GW; 
  REPORT(growth_a);
  REPORT(growth_re);
  REPORT(growth_repars);
  if(do_post_samp(5) == 1) ADREPORT(growth_re);

  // ---------------------------------------------------------------------
  // LAA re section
  Type nll_LAA = Type(0);

  if(LAA_re_model > 1) // random effects on LAA
  {
  
  Type sigma_LAA = exp(LAA_repars(0)); // first RE parameter
  Type rho_LAA_a = rho_trans(LAA_repars(1));
  Type rho_LAA_y = rho_trans(LAA_repars(2));  // second RE parameter
  Type Sigma_LAA;

	  if((LAA_re_model == 2) | (LAA_re_model == 5)) { // Only for ii_y and Ar1_y. CHECK THIS LATER, ONLY FOR WORK AR_y so far I think
		// likelihood of LAA deviations
		  Sigma_LAA = pow(pow(sigma_LAA,2) / ((1-pow(rho_LAA_y,2))*(1-pow(rho_LAA_a,2))),0.5);
		  nll_LAA += SCALE(SEPARABLE(AR1(rho_LAA_a),AR1(rho_LAA_y)), Sigma_LAA)(LAA_re); // must be array, not matrix!
		  SIMULATE if(simulate_state(5) == 1) if(sum(simulate_period) > 0) {
			array<Type> LAAre_tmp = LAA_re;
			SEPARABLE(AR1(rho_LAA_a),AR1(rho_LAA_y)).simulate(LAAre_tmp);
			LAAre_tmp = Sigma_LAA * LAAre_tmp;
			for(int y = 0; y < n_years_model + n_years_proj; y++){
			  if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
				for(int a = 0; a < n_ages; a++) LAA_re(y,a) = LAAre_tmp(y,a);
			  }
			}
		  }
	  } else { // for 2 and 3 options
        vector<Type> LAAre0 = LAA_re.matrix().row(0);
        Sigma_LAA = pow(pow(sigma_LAA,2) / (1-pow(rho_LAA_a,2)),0.5);
        nll_LAA += SCALE(AR1(rho_LAA_a), Sigma_LAA)(LAAre0);
        SIMULATE if(simulate_state(5) == 1) if(sum(simulate_period) > 0) {
          AR1(rho_LAA_a).simulate(LAAre0);
          for(int i = 0; i < LAAre0.size(); i++) LAAre0(i) = Sigma_LAA * LAAre0(i);
          for(int y = 0; y < n_years_model + n_years_proj; y++){
            for(int i = 0; i < LAAre0.size(); i++) LAA_re(y,i) = LAAre0(i);
          }
        }          
      }
	   
	if(do_post_samp.sum()==0){
		ADREPORT(sigma_LAA);
		ADREPORT(rho_LAA_a);
		ADREPORT(rho_LAA_y);
	}
		  
  }
  
  REPORT(nll_LAA);
  nll += nll_LAA; 
  REPORT(LAA_a);
  REPORT(LAA_re);
  REPORT(LAA_repars);
  if(do_post_samp(5) == 1) ADREPORT(LAA_re);

  // ---------------------------------------------------------------------
  // WAA re
  
  Type nll_WAA = Type(0);

	if(waa_type == 4) { // only if WAA deviates active
	
	  if(WAA_re_model > 1) // random effects on WAA
	  {
	  
	  Type sigma_WAA = exp(WAA_repars(0)); // first RE parameter
	  Type rho_WAA_a = rho_trans(WAA_repars(1));
	  Type rho_WAA_y = rho_trans(WAA_repars(2));  // second RE parameter
	  Type Sigma_WAA;

		  if((WAA_re_model == 2) | (WAA_re_model == 5)) { // 
			// likelihood of WAA deviations
			  Sigma_WAA = pow(pow(sigma_WAA,2) / ((1-pow(rho_WAA_y,2))*(1-pow(rho_WAA_a,2))),0.5);
			  nll_WAA += SCALE(SEPARABLE(AR1(rho_WAA_a),AR1(rho_WAA_y)), Sigma_WAA)(WAA_re); // must be array, not matrix!
			  SIMULATE if(simulate_state(5) == 1) if(sum(simulate_period) > 0) {
				array<Type> WAAre_tmp = WAA_re;
				SEPARABLE(AR1(rho_WAA_a),AR1(rho_WAA_y)).simulate(WAAre_tmp);
				WAAre_tmp = Sigma_WAA * WAAre_tmp;
				for(int y = 0; y < n_years_model + n_years_proj; y++){
				  if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
					for(int a = 0; a < n_ages; a++) WAA_re(y,a) = WAAre_tmp(y,a);
				  }
				}
			  }
		  } else { // for 2 and 3 options
			vector<Type> WAAre0 = WAA_re.matrix().row(0);
			Sigma_WAA = pow(pow(sigma_WAA,2) / (1-pow(rho_WAA_a,2)),0.5);
			nll_WAA += SCALE(AR1(rho_WAA_a), Sigma_WAA)(WAAre0);
			SIMULATE if(simulate_state(5) == 1) if(sum(simulate_period) > 0) {
			  AR1(rho_WAA_a).simulate(WAAre0);
			  for(int i = 0; i < WAAre0.size(); i++) WAAre0(i) = Sigma_WAA * WAAre0(i);
			  for(int y = 0; y < n_years_model + n_years_proj; y++){
				for(int i = 0; i < WAAre0.size(); i++) WAA_re(y,i) = WAAre0(i);
			  }
			}          
		  }
		   
		if(do_post_samp.sum()==0){
			ADREPORT(sigma_WAA);
			ADREPORT(rho_WAA_a);
			ADREPORT(rho_WAA_y);
		}
			  
	  }
	  
	} // waa_type
  
  REPORT(nll_WAA);
  nll += nll_WAA; 
  REPORT(WAA_a);
  REPORT(WAA_re);
  REPORT(WAA_repars);
  if(do_post_samp(5) == 1) ADREPORT(WAA_re);

  // ---------------------------------------------------------------------
  // Construct growth parameters per year
  array<Type> GW_par(n_years_model + n_years_proj,n_ages,n_growth_par); // array for growth parameters, either for 1 and 2
  for(int j = 0; j < n_growth_par; j++) { 
	  for(int y = 0; y < n_years_model; y++) { 
			for(int a = 0; a < n_ages; a++) { 
				GW_par(y,a,j) = exp(growth_a(j,0) + growth_re(y,a,j)); 
			}
	  }
   }

  // add to growth parameters in projection years
  if((do_proj == 1) & (growth_model == 1)){ 
  	for(int j = 0; j < n_growth_par; j++) { 
	  if(proj_GW_opt(j) == 2){
		  matrix<Type> GW_toavg(n_toavg,n_ages);
		  for(int a = 0; a < n_ages; a++){
			for(int i = 0; i < n_toavg; i++){
			  GW_toavg(i,a) = GW_par(avg_years_ind(i),a,j);
			}
		  }
		  vector<Type> GW_proj = GW_toavg.colwise().mean();
		  for(int y = n_years_model; y < n_years_model + n_years_proj; y++){
			for(int a = 0; a < n_ages; a++){
			  GW_par(y,a,j) = GW_proj(a);
			}
		  }
	  } else { // proj_GW_opt == 1
		for(int y = n_years_model; y < n_years_model + n_years_proj; y++) {
			for(int a = 0; a < n_ages; a++) { 		
				GW_par(y,a,j) = exp(growth_a(j,0) + growth_re(y,a,j)); 
			}
		}
	  }
	}
  }

  // Construct LAA parameters per year
  matrix<Type> LAA_par(n_years_model + n_years_proj,n_ages); 
  for(int y = 0; y < n_years_model; y++) { 
		for(int a = 0; a < n_ages; a++) { 
			LAA_par(y,a) = exp(LAA_a(a) + LAA_re(y,a)); 
		}
  }


  // add to LAA parameters in projection years
  if((do_proj == 1) & (growth_model == 2)){ 
	if(proj_GW_opt(0) == 2){
	  matrix<Type> GW_toavg(n_toavg,n_ages);
	  for(int a = 0; a < n_ages; a++){
		for(int i = 0; i < n_toavg; i++){
		  GW_toavg(i,a) = LAA_par(avg_years_ind(i),a);
		}
	  }
	  vector<Type> GW_proj = GW_toavg.colwise().mean();
	  for(int y = n_years_model; y < n_years_model + n_years_proj; y++){
		LAA_par.row(y) = GW_proj;
	  }
	} else { // proj_GW_opt == 1
		for(int y = n_years_model; y < n_years_model + n_years_proj; y++) {
			for(int a = 0; a < n_ages; a++) { 		
				LAA_par(y,a) = exp(LAA_a(a) + LAA_re(y,a)); 
			}
		}
	}
  }
  
  // add ecov effect on growth paramteres
  for(int j = 0; j < n_growth_par; j++) { 
	for(int i=0; i < n_Ecov; i++){
		if((Ecov_where(i,n_effects-2) == 1) & ((Ecov_where_subindex-1) == j)) {  // for growth
			for(int y = 0; y < n_years_model + n_years_proj; y++) {
				for(int a = 0; a < n_ages; a++) { 
					if(growth_model == 1) GW_par(y,a,j) *= exp(Ecov_lm(i,n_effects-2,y,0));
					if(growth_model == 2) LAA_par(y,a) *= exp(Ecov_lm(i,n_effects-2,y,a));
				}
			}
		}
    }
  }

  // ---------------------------------------------------------------------
  // Construct WAA parameters per year
  matrix<Type> WAA_par(n_years_model + n_years_proj,n_ages); 
  for(int y = 0; y < n_years_model; y++) { 
		for(int a = 0; a < n_ages; a++) { 
			WAA_par(y,a) = exp(WAA_a(a) + WAA_re(y,a)); 
		}
  }


  // add to WAA parameters in projection years
  if((do_proj == 1) & (waa_type == 4)){ 
	if(proj_WAA_opt(0) == 2){
	  matrix<Type> WAA_toavg(n_toavg,n_ages);
	  for(int a = 0; a < n_ages; a++){
		for(int i = 0; i < n_toavg; i++){
		  WAA_toavg(i,a) = WAA_par(avg_years_ind(i),a);
		}
	  }
	  vector<Type> WAA_proj = WAA_toavg.colwise().mean();
	  for(int y = n_years_model; y < n_years_model + n_years_proj; y++){
		WAA_par.row(y) = WAA_proj;
	  }
	} else { // proj_GW_opt == 1
		for(int y = n_years_model; y < n_years_model + n_years_proj; y++) {
			for(int a = 0; a < n_ages; a++) { 		
				WAA_par(y,a) = exp(WAA_a(a) + WAA_re(y,a)); 
			}
		}
	}
  }
  
  // Add Ecov for WAA



  // --------------------------------------------------------------------------
  // Calculate LW --  NEWG section
  Type nll_LW = Type(0);
  
  // For all growth parameters:
  vector<Type> sigma_LW(n_LW_par); // first RE parameter
  vector<Type> rho_LW_y(n_LW_par);  // second RE parameter

  for(int j = 0; j < n_LW_par; j++) {
  
  Type Sigma_LW = Type(0);

	  if((LW_re_model(j) == 2) | (LW_re_model(j) == 4)) { // Only for ii_y and Ar1_y. CHECK THIS LATER, ONLY FOR WORK AR_y so far I think
		
		sigma_LW(j) = exp(LW_repars(j,0));
        rho_LW_y(j) = rho_trans(LW_repars(j,1));  
		// likelihood of growth parameters deviations
			vector<Type> LWre0(n_years_model + n_years_proj);
			for(int y = 0; y < n_years_model + n_years_proj; y++) LWre0(y) = LW_re(y,0,j); 
			Sigma_LW = pow(pow(sigma_LW(j),2) / (1-pow(rho_LW_y(j),2)),0.5);
			nll_LW += SCALE(AR1(rho_LW_y(j)), Sigma_LW)(LWre0);
			SIMULATE if(simulate_state(6) == 1) if(sum(simulate_period) > 0) {
			  AR1(rho_LW_y(j)).simulate(LWre0);
			  for(int i = 0; i < LWre0.size(); i++) LWre0(i) = Sigma_LW * LWre0(i);
			  for(int y = 0; y < n_years_model + n_years_proj; y++){
				if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
				  for(int a = 0; a < n_ages; a++) LW_re(y,a,j) = LWre0(y);
				}
			  }
			}

		if(do_post_samp.sum()==0){
		  ADREPORT(sigma_LW);
		  ADREPORT(rho_LW_y);
		}
		
	  }
	  
	  if((LW_re_model(j) == 3) | (LW_re_model(j) == 5)) { // Only for ii_c and Ar1_c
		  
		sigma_LW(j) = exp(LW_repars(j,0));
        rho_LW_y(j) = rho_trans(LW_repars(j,1));  
		// likelihood of growth parameters deviations
			vector<Type> LWre0(n_years_model + n_years_proj + n_ages - 1);
			for(int i = 0; i < (n_ages - 1); i++) LWre0(i) = LW_re(0,n_ages - i - 1,j); // for cohorts at y = 0 except a = 0
			for(int i = (n_ages - 1); i < (n_years_model + n_years_proj + n_ages - 1); i++) LWre0(i) = LW_re(i-n_ages+1,0,j); // for cohorts y>=0 
			Sigma_LW = pow(pow(sigma_LW(j),2) / (1-pow(rho_LW_y(j),2)),0.5);
			nll_LW += SCALE(AR1(rho_LW_y(j)), Sigma_LW)(LWre0);
			SIMULATE if(simulate_state(6) == 1) if(sum(simulate_period) > 0) {
			  AR1(rho_LW_y(j)).simulate(LWre0);
			  for(int i = 0; i < LWre0.size(); i++) LWre0(i) = Sigma_LW * LWre0(i);
			  for(int y = 0; y < n_years_model + n_years_proj; y++){
				if(((simulate_period(0) == 1) & (y < n_years_model)) | ((simulate_period(1) == 1) & (y > n_years_model-1))){
				  for(int a = 0; a < n_ages; a++) {
					  if(y == 0) LW_re(y,n_ages-1-a,j) = LWre0(a); 
					  else LW_re(y,a,j) = LWre0(y-a+n_ages-1);
				  }
				}
			  }
			}

		if(do_post_samp.sum()==0){
		  ADREPORT(sigma_LW);
		  ADREPORT(rho_LW_y);
		}
		  
	  }
  
  }
  
  REPORT(nll_LW);
  nll += nll_LW; 
  REPORT(LW_a);
  REPORT(LW_re);
  REPORT(LW_repars);
  if(do_post_samp(6) == 1) ADREPORT(LW_re);

  // Construct LW parameters per year NEWG
  array<Type> LW_par(n_years_model + n_years_proj,n_ages,n_LW_par); // array for LW parameters
  for(int j = 0; j < n_LW_par; j++) { 
	  for(int y = 0; y < n_years_model; y++) { 
			for(int a = 0; a < n_ages; a++) { 
				LW_par(y,a,j) = exp(LW_a(j,0) + LW_re(y,a,j)); 
			}
	  }
   }


  // add to LW parameters in projection years
  if((do_proj == 1) & (waa_type > 1)){ 
  	for(int j = 0; j < n_LW_par; j++) { 
	  if(proj_LW_opt(j) == 2){
		  matrix<Type> LW_toavg(n_toavg,n_ages);
		  for(int a = 0; a < n_ages; a++){
			for(int i = 0; i < n_toavg; i++){
			  LW_toavg(i,a) = LW_par(avg_years_ind(i),a,j);
			}
		  }
		  vector<Type> LW_proj = LW_toavg.colwise().mean();
		  for(int y = n_years_model; y < n_years_model + n_years_proj; y++){
			for(int a = 0; a < n_ages; a++){
			  LW_par(y,a,j) = LW_proj(a);
			}
		  }
	  } else { // proj_GW_opt == 1
		for(int y = n_years_model; y < n_years_model + n_years_proj; y++) {
			for(int a = 0; a < n_ages; a++) { 		
				LW_par(y,a,j) = exp(LW_a(j,0) + LW_re(y,a,j)); 
			}
		}
	  }
	}
  }
  
  // add ecov effect on LW paramteres
  for(int j = 0; j < n_LW_par; j++) { 
	for(int i=0; i < n_Ecov; i++){
		if((Ecov_where(i,n_effects-1) == 1) & ((Ecov_where_subindex-1) == j)) {  // for LW
			for(int y = 0; y < n_years_model + n_years_proj; y++) {
				for(int a = 0; a < n_ages; a++) { 
					LW_par(y,a,j) *= exp(Ecov_lm(i,n_effects-1,y,0));
				}
			}
		}
    }
  }  
  
  // --------------------------------------------------------------------------
  // Calculate mean-LAA, SDAA, and transition matrix, for all years: NEWG
  // Some information I need for transition matrix: 
  Type len_bin = lengths(1) - lengths(0); // input should have standardized length bin
  Type Lminp = min(lengths) + len_bin*0.5;
  Type Lmaxp = max(lengths) - len_bin*0.5;
  Type Fac1 = 0.0;
  Type Fac2 = 0.0;
  Type Ll1p = 0.0;
  Type Llp = 0.0;
  for(int y = 0; y < n_years_model + n_years_proj; y++)
  {
	  for(int a = 0; a < n_ages; a++)
	  {
			if(growth_model == 1) {
				if(y == 0) {
					LAA(y,a) = GW_par(y,a,1) + (GW_par(y,a,2) - GW_par(y,a,1)) * exp(-GW_par(y,a,0)*a); // for growth_model = 1
				} else {
					if(a == 0) { 
						LAA(y,a) = GW_par(y,a,2);  
					} else {
						LAA(y,a) = LAA(y-1,a-1) + (LAA(y-1,a-1) - GW_par(y-1,a-1,1))*(exp(-GW_par(y-1,a-1,0)) - 1.0); // use growth parameters y-1 and a-1 because it is jan1
					}
				}
			}
			if(growth_model == 2) LAA(y,a) = LAA_par(y,a); // for growth_model = 2
			
			if(growth_model == 1) SDAA(y,a) = ( GW_par(y,a,3) + ((GW_par(y,a,4) - GW_par(y,a,3))/(n_ages - 1.0))*a )*LAA(y,a);  //
			if(growth_model == 2) SDAA(y,a) = ( GW_par(y,a,0) + ((GW_par(y,a,1) - GW_par(y,a,0))/(n_ages - 1.0))*a )*LAA(y,a);  //
			
			// pred_LAA(0,a) = LAA(0,a); // predicted mean length at age, is it necessary?
			for(int l = 0; l < n_lengths; l++) {
				
				if(l == 0) { 
					Fac1 = (Lminp - LAA(y,a))/SDAA(y,a);
					phi_mat(y,l,a) = pnorm(Fac1);  
				} else {
					if(l == (n_lengths-1)) { 
						Fac1 = (Lmaxp - LAA(y,a))/SDAA(y,a);
						phi_mat(y,l,a) = 1.0 - pnorm(Fac1);  
					} else { 
						Ll1p = lengths(l) + len_bin*0.5;
						Llp = lengths(l) - len_bin*0.5;
						Fac1 = (Ll1p - LAA(y,a))/SDAA(y,a);
						Fac2 = (Llp - LAA(y,a))/SDAA(y,a);
						phi_mat(y,l,a) = pnorm(Fac1) - pnorm(Fac2);  
					}
				}
				
			}
	  }
  }

  // --------------------------------------------------------------------------
  // Weight at age calculations:
  	Type sum_wt = 0;
	Type sum_wt_fleet = 0;
	Type sum_wt_index = 0;
	matrix<Type> watl(n_years_model + n_years_proj, n_lengths);
	array<Type> nll_waa(waa.dim(0), n_years_model, n_ages);
	nll_waa.setZero();
  if(waa_type == 1) {
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
	  if(waa_type == 2) { // waa not provided, so use LW parameters, waa array will be replaced by pred_waa
		  for(int y = 0; y < n_years_model + n_years_proj; y++) {
			int yuse = y;
			int y_1 = y + 1;
			if(y > n_years_model - 1) yuse = n_years_model -1; //some things only go up to n_years_model-1
			if(y == (n_years_model + n_years_proj - 1)) y_1 = y;
			
			// For Jan-1
			for(int a = 0; a < n_ages; a++) { 
				sum_wt = 0;
				for(int l = 0; l < n_lengths; l++) {
					watl(y,l) = LW_par(y,a,0)*pow(lengths(l), LW_par(y,a,1));
					if(phi_matrix_info == 0) sum_wt += phi_mat(y,l,a)*watl(y,l);
					else sum_wt += phi_matrix_input(waa_pointer_jan1-1,l,a)*watl(y,l); 
				}
				pred_waa(waa_pointer_jan1 - 1,y,a) = sum_wt; // jan-1st
			}
			
			// For SSB
			fracyr_phi_mat = pred_LAA(vector<Type>(LAA.row(y)), vector<Type>(SDAA.row(y)), vector<Type>(LAA.row(y_1)), GW_par, lengths, y, fracyr_SSB(yuse), growth_model);
			for(int a = 0; a < n_ages; a++) { 
				sum_wt_ssb = 0;
				for(int l = 0; l < n_lengths; l++) {
					if(phi_matrix_info == 0) sum_wt_ssb += fracyr_phi_mat(l,a)*watl(y,l);
					else sum_wt_ssb += phi_matrix_input(waa_pointer_ssb-1,l,a)*watl(y,l); 
				}
				pred_waa(waa_pointer_ssb - 1,y,a) = sum_wt_ssb; // SSB
			}

			// For fleets
			for(int f = 0; f < n_fleets; f++) {
				fracyr_phi_mat = pred_LAA(vector<Type>(LAA.row(y)), vector<Type>(SDAA.row(y)), vector<Type>(LAA.row(y_1)), GW_par, lengths, y, fracfleet(yuse,f), growth_model);
				for(int a = 0; a < n_ages; a++) { 
					sum_wt_fleet = 0;
					for(int l = 0; l < n_lengths; l++) {
						if(phi_matrix_info == 0) sum_wt_fleet += fracyr_phi_mat(l,a)*watl(y,l);
						else sum_wt_fleet += phi_matrix_input(waa_pointer_fleets(f)-1,l,a)*watl(y,l);
					}
					pred_waa(waa_pointer_fleets(f)-1,y,a) = sum_wt_fleet; 
					pred_waa(waa_pointer_totcatch-1,y,a) = sum_wt_fleet; // for total catch
				}
			}
			
			// For indices
			for(int i = 0; i < n_indices; i++) {
				fracyr_phi_mat = pred_LAA(vector<Type>(LAA.row(y)), vector<Type>(SDAA.row(y)), vector<Type>(LAA.row(y_1)), GW_par, lengths, y, fracyr_indices(yuse,i), growth_model);
				for(int a = 0; a < n_ages; a++) { 
					sum_wt_index = 0;
					for(int l = 0; l < n_lengths; l++) {
						if(phi_matrix_info == 0) sum_wt_index += fracyr_phi_mat(l,a)*watl(y,l);
						else sum_wt_index += phi_matrix_input(waa_pointer_indices(i)-1,l,a)*watl(y,l);
					}
					pred_waa(waa_pointer_indices(i)-1,y,a) = sum_wt_index; // for indices	
				}
			}
			
		  }
	  } else { // waa provided and LW parameters used, likelihood function added
		if(waa_type == 3) {
			  for(int y = 0; y < n_years_model + n_years_proj; y++) { // 
				int yuse = y;
				int y_1 = y + 1;
				if(y > n_years_model - 1) yuse = n_years_model -1; //some things only go up to n_years_model-1
				if(y == (n_years_model + n_years_proj - 1)) y_1 = y;
				
				// For Jan-1
				for(int a = 0; a < n_ages; a++) { 
					sum_wt = 0;
					for(int l = 0; l < n_lengths; l++) {
						watl(y,l) = LW_par(y,a,0)*pow(lengths(l), LW_par(y,a,1));
						if(phi_matrix_info == 0) sum_wt += phi_mat(y,l,a)*watl(y,l);
						else sum_wt += phi_matrix_input(waa_pointer_jan1-1,l,a)*watl(y,l); 
					}
					pred_waa(waa_pointer_jan1 - 1,y,a) = sum_wt; // jan-1st
				}
				
				// For SSB
				fracyr_phi_mat = pred_LAA(vector<Type>(LAA.row(y)), vector<Type>(SDAA.row(y)), vector<Type>(LAA.row(y_1)), GW_par, lengths, y, fracyr_SSB(yuse), growth_model);
				for(int a = 0; a < n_ages; a++) { 
					sum_wt_ssb = 0;
					for(int l = 0; l < n_lengths; l++) {
						if(phi_matrix_info == 0) sum_wt_ssb += fracyr_phi_mat(l,a)*watl(y,l);
						else sum_wt_ssb += phi_matrix_input(waa_pointer_ssb-1,l,a)*watl(y,l); 
					}
					pred_waa(waa_pointer_ssb - 1,y,a) = sum_wt_ssb; // SSB
				}

				// For fleets
				for(int f = 0; f < n_fleets; f++) {
					fracyr_phi_mat = pred_LAA(vector<Type>(LAA.row(y)), vector<Type>(SDAA.row(y)), vector<Type>(LAA.row(y_1)), GW_par, lengths, y, fracfleet(yuse,f), growth_model);
					for(int a = 0; a < n_ages; a++) { 
						sum_wt_fleet = 0;
						for(int l = 0; l < n_lengths; l++) {
							if(phi_matrix_info == 0) sum_wt_fleet += fracyr_phi_mat(l,a)*watl(y,l);
							else sum_wt_fleet += phi_matrix_input(waa_pointer_fleets(f)-1,l,a)*watl(y,l);
						}
						pred_waa(waa_pointer_fleets(f)-1,y,a) = sum_wt_fleet; 
						pred_waa(waa_pointer_totcatch-1,y,a) = sum_wt_fleet; // for total catch
						if((y < n_years_model) & (waa_cv(waa_pointer_fleets(f) - 1,y,a) > 0) & (use_catch_waa(y,f) == 1)) { // add here for totalcatch if required
							nll_waa(waa_pointer_fleets(f) - 1,y,a) += get_waa_ll(waa(waa_pointer_fleets(f) - 1,y,a), pred_waa(waa_pointer_fleets(f) - 1,y,a), waa_cv(waa_pointer_fleets(f) - 1,y,a)); 
						}
					}
				}
				
				// For indices
				for(int i = 0; i < n_indices; i++) {
					fracyr_phi_mat = pred_LAA(vector<Type>(LAA.row(y)), vector<Type>(SDAA.row(y)), vector<Type>(LAA.row(y_1)), GW_par, lengths, y, fracyr_indices(yuse,i), growth_model);
					for(int a = 0; a < n_ages; a++) { 
						sum_wt_index = 0;
						for(int l = 0; l < n_lengths; l++) {
							if(phi_matrix_info == 0) sum_wt_index += fracyr_phi_mat(l,a)*watl(y,l);
							else sum_wt_index += phi_matrix_input(waa_pointer_indices(i)-1,l,a)*watl(y,l);
						}
						pred_waa(waa_pointer_indices(i)-1,y,a) = sum_wt_index; // for indices	
						if((y < n_years_model) & (waa_cv(waa_pointer_indices(i) - 1,y,a) > 0) & (use_index_waa(y,i) == 1)) {
							nll_waa(waa_pointer_indices(i) - 1,y,a) += get_waa_ll(waa(waa_pointer_indices(i) - 1,y,a), pred_waa(waa_pointer_indices(i) - 1,y,a), waa_cv(waa_pointer_indices(i) - 1,y,a)); 
						}
					}
				}
				
			  }  // loop only for year_model
			  nll += nll_waa.sum();	
		} else {
			// waa_type == 4
			for(int y = 0; y < n_years_model + n_years_proj; y++) {
				int yuse = y;
				int y_1 = y + 1;
				if(y > n_years_model - 1) yuse = n_years_model -1; //some things only go up to n_years_model-1
				if(y == (n_years_model + n_years_proj - 1)) y_1 = y;
				
				// For Jan-1
				for(int a = 0; a < n_ages; a++) { 
					pred_waa(waa_pointer_jan1 - 1,y,a) = WAA_par(y,a); // jan-1st
				}
				// For SSB
				fracyr_WAA = get_fracyr_WAA(vector<Type>(WAA_par.row(y)), vector<Type>(WAA_par.row(y_1)), fracyr_SSB(yuse));
				for(int a = 0; a < n_ages; a++) { 
					pred_waa(waa_pointer_ssb - 1,y,a) = fracyr_WAA(a); // SSB
				}
				// For fleets
				for(int f = 0; f < n_fleets; f++) {
					fracyr_WAA = get_fracyr_WAA(vector<Type>(WAA_par.row(y)), vector<Type>(WAA_par.row(y_1)), fracfleet(yuse,f));
					for(int a = 0; a < n_ages; a++) { 
						pred_waa(waa_pointer_fleets(f)-1,y,a) = fracyr_WAA(a); 
						pred_waa(waa_pointer_totcatch-1,y,a) = fracyr_WAA(a); // for total catch
						if((y < n_years_model) & (waa_cv(waa_pointer_fleets(f) - 1,y,a) > 0) & (use_catch_waa(y,f) == 1)) { // add here for totalcatch if required
							nll_waa(waa_pointer_fleets(f) - 1,y,a) += get_waa_ll(waa(waa_pointer_fleets(f) - 1,y,a), pred_waa(waa_pointer_fleets(f) - 1,y,a), waa_cv(waa_pointer_fleets(f) - 1,y,a)); 
						}
					}
				}
				// For indices
				for(int i = 0; i < n_indices; i++) {
					fracyr_WAA = get_fracyr_WAA(vector<Type>(WAA_par.row(y)), vector<Type>(WAA_par.row(y_1)), fracyr_indices(yuse,i));
					for(int a = 0; a < n_ages; a++) { 
						pred_waa(waa_pointer_indices(i)-1,y,a) = fracyr_WAA(a); // for indices	
						if((y < n_years_model) & (waa_cv(waa_pointer_indices(i) - 1,y,a) > 0) & (use_index_waa(y,i) == 1)) {
							nll_waa(waa_pointer_indices(i) - 1,y,a) += get_waa_ll(waa(waa_pointer_indices(i) - 1,y,a), pred_waa(waa_pointer_indices(i) - 1,y,a), waa_cv(waa_pointer_indices(i) - 1,y,a)); 
						}
					}
				}				

			}
			
		} 
	  }
  }
  REPORT(pred_waa); // print predicted waa matrix	  
  REPORT(nll_waa);

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
  
  // Transform selex-at-len to selex-at-age:
  for(int b = 0; b < n_selblocks; b++) {
	  if(selblock_models(b) < 6) { // for age models
		  selAA(b) = selAL(b);  // same as calculated in selAL
	  } else { // transform for all years
	      matrix<Type> matemp(n_years_model, n_ages);
		  for(int y = 0; y < n_years_model; y++){
			  for(int a = 0; a < n_ages; a++) {
			  Type sumSelex = 0.0;
				  for(int l = 0; l < n_lengths; l++) {
					if(phi_matrix_info == 0) sumSelex += phi_mat(y,l,a)*selAL(b)(y,l);
					else sumSelex += phi_matrix_input(waa_pointer_jan1-1,l,a)*selAL(b)(y,l);
				  }
			  matemp(y,a) = sumSelex;
			  }
		  }
		  selAA(b) = matemp;
	  }
  }
  
  // Construct survey catchability-at-age (QAA)
  for(int i = 0; i < n_indices; i++)
  {
  // add ecov effect on M (by year, shared across ages)
    for(int y = 0; y < n_years_model; y++)
    {
      for(int a = 0; a < n_ages; a++) QAA(y,i,a) = q(y,i) * selAA(selblock_pointer_indices(y,i)-1)(y,a);
    }
    //just use last years selectivity for now
    if(do_proj == 1) for(int y = n_years_model; y < n_years_model + n_years_proj; y++) for(int a = 0; a < n_ages; a++) {
      QAA(y,i,a) = q(y,i) * selAA(selblock_pointer_indices(n_years_model-1,i)-1)(n_years_model-1,a);
    }
  }
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
  REPORT(QAA);

  // Construct fishing mortality-at-age (FAA)
  FAA_tot.setZero();
  for(int f = 0; f < n_fleets; f++)
  {
    log_F(0,f) = log_F1(f);
    F(0,f) = exp(log_F(0,f));
    for(int a = 0; a < n_ages; a++)
    {
      FAA(0,f,a) = F(0,f) * selAA(selblock_pointer_fleets(0,f)-1)(0,a);
      log_FAA(0,f,a) = log(FAA(0,f,a));
      FAA_tot(0,a) = FAA_tot(0,a) + FAA(0,f,a);
    }
    for(int y = 1; y < n_years_model; y++)
    {
      log_F(y,f) = log_F(y-1,f) + F_devs(y-1,f);
      F(y,f) = exp(log_F(y,f));
      for(int a = 0; a < n_ages; a++)
      {
        FAA(y,f,a) = F(y,f) * selAA(selblock_pointer_fleets(y,f)-1)(y,a);
        log_FAA(y,f,a) = log(FAA(y,f,a));
        FAA_tot(y,a) = FAA_tot(y,a) + FAA(y,f,a);
      }
    }
  }
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
    SSB(0) += NAA(0,a) * pred_waa(waa_pointer_ssb-1,0,a) * mature(0,a) * exp(-ZAA(0,a)*fracyr_SSB(0));
    pred_NAA(0,a) = NAA(0,a);
  }

  // get SPR0
  vector<Type> M(n_ages), sel(n_ages), mat(n_ages), waacatch(n_ages), waassb(n_ages), log_SPR0(n_years_model + n_years_proj);
  int na = n_years_model + n_years_proj;
  vector<Type> log_SR_a(na), log_SR_b(na), SR_h(na), SR_R0(na);
  for(int y = 0; y < n_years_model + n_years_proj; y++)
  {
    for(int a = 0; a < n_ages; a++)
    {
      M(a) = MAA(y,a);
      waassb(a) = pred_waa(waa_pointer_ssb-1,y,a);
      mat(a) = mature(y,a);
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
		waacatch = get_waa_y(pred_waa, y, n_ages, waa_pointer_totcatch);
		waassb = get_waa_y(pred_waa, y, n_ages, waa_pointer_ssb);
        //n_fleets x n_ages: projected full F is sum of (means across years at age) across fleets 
        matrix<Type> FAA_proj = get_F_proj(y, n_fleets, proj_F_opt, FAA, NAA, MAA, mature, waacatch, waassb, fracyr_SSB, 
          log_SPR0, avg_years_ind, n_years_model, which_F_age, percentSPR, proj_Fcatch, percentFXSPR, F_proj_init(y-n_years_model), 
          log_SR_a, log_SR_b, recruit_model, percentFMSY);
        FAA_tot.row(y) = FAA_proj.colwise().sum();
        for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) FAA(y,f,a) = FAA_proj(f,a);
        //FAA_tot.row(y) = get_F_proj(y, proj_F_opt, FAA_tot, NAA, MAA, mature, waacatch, waassb, fracyr_SSB, log_SPR0, avg_years_ind, n_years_model,
        // which_F_age, percentSPR, proj_Fcatch);
        ZAA.row(y) = FAA_tot.row(y) + MAA.row(y);
      }
    } // end proj F
    SSB(y) = get_SSB(NAA,ZAA,pred_waa, mature,y, waa_pointer_ssb, fracyr_SSB);
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
    matrix<Type> sims = sim_pop(NAA_devs, recruit_model, mean_rec_pars, SSB, NAA, log_SR_a, log_SR_b, Ecov_where, Ecov_how, Ecov_lm, 
      n_NAA_sigma, do_proj, proj_F_opt, FAA, FAA_tot, MAA, mature, pred_waa, waa_pointer_totcatch, waa_pointer_ssb, fracyr_SSB, log_SPR0, 
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

  // Catch data likelihood
  matrix<Type> nll_agg_catch(n_years_model,n_fleets), nll_catch_acomp(n_years_model,n_fleets), nll_catch_lcomp(n_years_model,n_fleets), agg_catch_proj(n_years_proj,n_fleets);
  array<Type> nll_catch_caal(n_years_model,n_fleets,n_lengths);
  array<Type> catch_paa_proj(n_fleets, n_years_proj, n_ages);
  array<Type> catch_pal_proj(n_fleets, n_years_proj, n_lengths); // NEWG
  array<Type> catch_caal_proj(n_fleets,n_years_proj,n_lengths, n_ages); // NEWG
  nll_agg_catch.setZero();
  nll_catch_acomp.setZero();
  nll_catch_lcomp.setZero(); // NEWG
  nll_catch_caal.setZero();
  vector<Type> lsum(n_lengths);
  for(int y = 0; y < n_years_model+n_years_proj; y++)
  {
    //for now just use uncertainty from last year of catch
    int usey = y;
	int y_1 = y + 1;
    if(y > n_years_model-1) usey = n_years_model-1;
	if(y == (n_years_model + n_years_proj - 1)) y_1 = y;
    //int acomp_par_count = 0;
	fracyr_phi_mat = pred_LAA(vector<Type>(LAA.row(y)), vector<Type>(SDAA.row(y)), vector<Type>(LAA.row(y_1)), GW_par, lengths, y, fracfleet, growth_model); // only works for growth_model = 1 so far
	for(int f = 0; f < n_fleets; f++)
    {
	  lsum.setZero();
      pred_catch(y,f) = 0.0;
      Type tsum = 0.0;
      for(int a = 0; a < n_ages; a++)
      {
        pred_CAA(y,f,a) =  NAA(y,a) * FAA(y,f,a) * (1 - exp(-ZAA(y,a)))/ZAA(y,a);
        pred_catch(y,f) += pred_waa(waa_pointer_fleets(f)-1,y,a) * pred_CAA(y,f,a);
        tsum += pred_CAA(y,f,a);
		// NEWG: is there a more efficient way to do this?:
		for(int l = 0; l < n_lengths; l++) {
			if(phi_matrix_info == 0) pred_CAAL(y,f,l,a) = pred_CAA(y,f,a) * fracyr_phi_mat(l,a);
			else pred_CAAL(y,f,l,a) = pred_CAA(y,f,a) * phi_matrix_input(waa_pointer_fleets(f)-1,l,a);
			lsum(l) += pred_CAAL(y,f,l,a);
		}
      }
	  for(int l = 0; l < n_lengths; l++) pred_CAL(y,f,l) = lsum(l); // NEWG
		       
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
        for(int a = 0; a < n_ages; a++){
          pred_catch_paa(y,f,a) = pred_CAA(y,f,a)/tsum;
          t_pred_paa(a) = pred_catch_paa(y,f,a);
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
          vector<int> lens_obs_y = lensvec.segment(keep_Cpal(f,y,0), keep_Cpal(f,y,1));
          nll_catch_lcomp(y,f) -= get_acomp_ll(tf_pal_obs, t_pred_pal, catch_NeffL(y,f), lens_obs_y, len_comp_model_fleets(f), 
            vector<Type>(catch_pal_pars.row(f)), keep.segment(keep_Cpal(f,y,0),keep_Cpal(f,y,1)), do_osa, pal_obs_y);
        }
        SIMULATE if(simulate_data(0) == 1) if(use_catch_pal(usey,f) == 1){
          if((simulate_period(0) == 1) & (y < n_years_model)) //model years
          {
            for(int l = 0; l < n_lengths; l++) pal_obs_y(l) = catch_pal(f,y,l);
            vector<int> lens_obs_y = lensvec.segment(keep_Cpal(f,y,0), keep_Cpal(f,y,1));
            vector<Type> tf_pal_obs = sim_acomp(t_pred_pal, catch_NeffL(y,f), lens_obs_y, len_comp_model_fleets(f), vector<Type>(catch_pal_pars.row(f)));
            obsvec.segment(keep_Cpal(f,y,0),keep_Cpal(f,y,1)) = tf_pal_obs;
            pal_obs_y = make_paa(tf_pal_obs, len_comp_model_fleets(f), lens_obs_y, pal_obs_y);
            for(int l = 0; l < n_lengths; l++) catch_pal(f,y,l) = pal_obs_y(l);
          }
          if((simulate_period(1) == 1) & (y > n_years_model - 1)) //projection years
          {
            for(int l = 0; l < n_lengths; l++) pal_obs_y(l) = 1/n_lengths; //only needed for LN obs to tell where the last non-zero age class is. No zeros in projections.
            vector<int> lens_obs_y(n_lengths);
            for(int l = 0; l < n_lengths; l++) lens_obs_y(l) = l;
            vector<Type> tf_pal_obs = sim_acomp(t_pred_pal, catch_NeffL(usey,f), lens_obs_y, len_comp_model_fleets(f), vector<Type>(catch_pal_pars.row(f)));
            pal_obs_y = make_paa(tf_pal_obs, len_comp_model_fleets(f), lens_obs_y, pal_obs_y);
            for(int l = 0; l < n_lengths; l++) catch_pal_proj(f,y-n_years_model,l) = pal_obs_y(l);
          }
        }
      }
	  
	  // CAAL nll 
	  if(any_fleet_caal(f) == 1){
		  for(int l = 0; l < n_lengths; l++) {
			vector<Type> paa_obs_y(n_ages);
			paa_obs_y.setZero();
			for(int a = 0; a < n_ages; a++){
			  pred_catch_paa(y,f,l,a) = pred_CAAL(y,f,l,a)/lsum(l);
			  t_pred_paa(a) = pred_catch_caal(y,f,l,a);
			}
			if(y < n_years_model) if(use_catch_caal(y,f) == 1) {
			  for(int a = 0; a < n_ages; a++) paa_obs_y(a) = catch_caal(f,y,l,a);
			  //NB: indexing in obsvec MUST be: keep_Cpaa(i,y,0),...,keep_Cpaa(i,y,0) + keep_Cpaa(i,y,1) - 1
			  //keep_Cpaa(i,y,0) is first val, keep_Cpaa(i,y,1) is the length of the vector
			  vector<Type> tf_paa_obs = obsvec.segment(keep_Ccaal(f,y,l,0), keep_Ccaal(f,y,l,1));
			  vector<int> ages_obs_y = agesvec.segment(keep_Ccaal(f,y,l,0), keep_Ccaal(f,y,l,1));
			  nll_catch_caal(y,f,l) -= get_acomp_ll(tf_paa_obs, t_pred_paa, catch_caal_Neff(y,f,l), ages_obs_y, age_comp_model_fleets(f), 
				vector<Type>(catch_paa_pars.row(f)), keep.segment(keep_Ccaal(f,y,l,0),keep_Ccaal(f,y,l,1)), do_osa, paa_obs_y);
			}
			SIMULATE if(simulate_data(0) == 1) if(use_catch_caal(usey,f) == 1){
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
		  }
      }
	 
    }
  }
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(0) == 1) REPORT(agg_catch);
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(0) == 1) REPORT(catch_paa);
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(0) == 1) REPORT(catch_pal); // NEWG
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(0) == 1) REPORT(catch_caal); // NEWG
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(1) == 1) REPORT(agg_catch_proj);
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(1) == 1) REPORT(catch_paa_proj);
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(1) == 1) REPORT(catch_pal_proj); // NEWG
  SIMULATE if(simulate_data(0) == 1) if(simulate_period(1) == 1) REPORT(catch_caal_proj); // NEWG
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
  array<Type> index_pal_proj(n_indices,n_years_proj,n_lengths); // NEWG
  array<Type> index_caal_proj(n_indices,n_years_proj,n_lengths, n_ages); // NEWG
  nll_agg_indices.setZero();
  nll_index_acomp.setZero();
  nll_index_lcomp.setZero();
  nll_index_caal.setZero();
  pred_indices.setZero();
  vector<Type> lsumI(n_lengths);
  for(int y = 0; y < n_years_model + n_years_proj; y++)
  {
    int usey = y;
	int y_1 = y + 1;
    if(y > n_years_model - 1) usey = n_years_model -1; //some things only go up to n_years_model-1
	if(y == (n_years_model + n_years_proj - 1)) y_1 = y;
    //int acomp_par_count = 0;
    for(int i = 0; i < n_indices; i++)
    {
	  lsumI.setZero();
      Type tsum = 0.0;
	  fracyr_phi_mat = pred_LAA(vector<Type>(LAA.row(y)), vector<Type>(SDAA.row(y)), vector<Type>(LAA.row(y_1)), GW_par, lengths, y, fracyr_indices(usey,i), growth_model); // only works for growth_model = 1 so far
      for(int a = 0; a < n_ages; a++)
      {
        pred_IAA(y,i,a) =  NAA(y,a) * QAA(y,i,a) * exp(-ZAA(y,a) * fracyr_indices(usey,i));
        if(units_indices(i) == 1) pred_indices(y,i) += pred_waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        else pred_indices(y,i) += pred_IAA(y,i,a);
		
		// NEWG: is there a more efficient way to do this?:
		for(int l = 0; l < n_lengths; l++) {
			if(phi_matrix_info == 0) pred_IAAL(y,i,l,a) = pred_IAA(y,i,a) * fracyr_phi_mat(l,a); // Only numbers allowed so far
			else pred_IAAL(y,i,l,a) = pred_IAA(y,i,a) * phi_matrix_input(waa_pointer_indices(i)-1,l,a);
			lsumI(l) += pred_IAAL(y,i,l,a);
		}
      }
	  for(int l = 0; l < n_lengths; l++) pred_IAL(y,i,l) = lsumI(l); // NEWG
      
	  for(int a = 0; a < n_ages; a++)
      {
        if(units_index_paa(i) == 1) pred_IAA(y,i,a) = pred_waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        tsum += pred_IAA(y,i,a);
      }

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
        for(int a = 0; a < n_ages; a++)
        {
          pred_index_paa(y,i,a) = pred_IAA(y,i,a)/tsum;
          t_pred_paa(a) = pred_index_paa(y,i,a);
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
          pred_index_pal(y,i,l) = pred_IAL(y,i,l)/lsumI.sum();
          t_pred_pal(l) = pred_index_pal(y,i,l);
        }
        if(y < n_years_model) if(use_index_pal(y,i) == 1) {
          for(int l = 0; l < n_lengths; l++) pal_obs_y(l) = index_pal(i,y,l);
          //NB: indexing in obsvec MUST be: keep_Cpaa(i,y,0),...,keep_Cpaa(i,y,0) + keep_Cpaa(i,y,1) - 1
          //keep_Cpaa(i,y,0) is first val, keep_Cpaa(i,y,1) is the length of the vector
          vector<Type> tf_pal_obs = obsvec.segment(keep_Ipal(i,y,0), keep_Ipal(i,y,1));
          vector<int> lens_obs_y = lensvec.segment(keep_Ipal(i,y,0), keep_Ipal(i,y,1));
          nll_index_lcomp(y,i) -= get_acomp_ll(tf_pal_obs, t_pred_pal, index_NeffL(y,i), lens_obs_y, len_comp_model_indices(i), 
            vector<Type>(index_pal_pars.row(i)), keep.segment(keep_Ipal(i,y,0),keep_Ipal(i,y,1)), do_osa, pal_obs_y);
        }
        SIMULATE if(simulate_data(1) == 1) if(use_index_pal(usey,i) == 1){
          if((simulate_period(0) == 1) & (y < n_years_model)) //model years
          {
            for(int l = 0; l < n_lengths; l++) pal_obs_y(l) = index_pal(i,y,l);
            vector<int> lens_obs_y = lensvec.segment(keep_Ipal(i,y,0), keep_Ipal(i,y,1));
            vector<Type> tf_pal_obs = sim_acomp(t_pred_pal, index_NeffL(usey,i), lens_obs_y, len_comp_model_indices(i), vector<Type>(index_pal_pars.row(i)));
            obsvec.segment(keep_Ipal(i,y,0),keep_Ipal(i,y,1)) = tf_pal_obs;
            pal_obs_y = make_paa(tf_pal_obs, len_comp_model_indices(i), lens_obs_y, pal_obs_y);
            for(int l = 0; l < n_lengths; l++) index_pal(i,y,l) = pal_obs_y(l);
          }
          if((simulate_period(1) == 1) & (y > n_years_model - 1)) //projection years
          {
            for(int l = 0; l < n_lengths; l++) pal_obs_y(l) = 1/n_lengths; //only needed for LN obs to tell where the last non-zero age class is. No zeros in projections.
            vector<int> lens_obs_y(n_lengths);
		  for(int l = 0; l < n_lengths; l++) lens_obs_y(l) = l;
            vector<Type> tf_pal_obs = sim_acomp(t_pred_pal, index_NeffL(usey,i), lens_obs_y, len_comp_model_indices(i), vector<Type>(index_pal_pars.row(i)));
            pal_obs_y = make_paa(tf_pal_obs, len_comp_model_indices(i), lens_obs_y, pal_obs_y);
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
			for(int a = 0; a < n_ages; a++)
			{
				pred_index_caal(y,i,l,a) = pred_IAAL(y,i,l,a)/lsumI(l);
				t_pred_paa(a) = pred_index_caal(y,i,l,a); 
			}
			if(y < n_years_model) if(use_index_caal(y,i) == 1) {
			  for(int a = 0; a < n_ages; a++) paa_obs_y(a) = index_caal(i,y,l,a);
			  //NB: indexing in obsvec MUST be: keep_Ipaa(i,y,0),...,keep_Ipaa(i,y,0) + keep_Ipaa(i,y,1) - 1
			  //keep_Ipaa(i,y,0) is first val, keep_Ipaa(i,y,1) is the length of the vector
			  vector<Type> tf_paa_obs = obsvec.segment(keep_Icaal(i,y,l,0), keep_Icaal(i,y,l,1));
			  vector<int> ages_obs_y = agesvec.segment(keep_Icaal(i,y,l,0), keep_Icaal(i,y,l,1));
			  nll_index_caal(y,i,l) -= get_acomp_ll(tf_paa_obs, t_pred_paa, index_caal_Neff(y,i,l), ages_obs_y, age_comp_model_indices(i), 
				vector<Type>(index_paa_pars.row(i)), keep.segment(keep_Icaal(i,y,l,0),keep_Icaal(i,y,l,1)), do_osa, paa_obs_y);
			}
			SIMULATE if(simulate_data(1) == 1) if(use_index_caal(usey,i) == 1){
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
    vector<Type> catch_proj(n_years_proj), log_catch_proj(n_years_proj);
    matrix<Type> CAA_proj(n_years_proj, n_ages), catch_fleet_proj(n_years_proj, n_fleets);
    array<Type> CAA_fleet_proj(n_fleets, n_years_proj, n_ages);
    catch_proj.setZero();
    for(int i = 0; i < n_years_proj; i++){
      int yi = i + n_years_model;
      for(int a = 0; a < n_ages; a++){
        CAA_proj(i,a) =  NAA(yi,a) * FAA_tot(yi,a) * (1 - exp(-ZAA(yi,a)))/ZAA(yi,a);
        waacatch(a) = pred_waa(waa_pointer_totcatch-1, yi, a);
        catch_proj(i) += waacatch(a) * CAA_proj(i,a);
      }
      log_catch_proj(i) = log(catch_proj(i) + Type(1.0e-15));
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
  matrix<Type> SPR_res = get_SPR_res(MAA, FAA_tot, which_F_age, pred_waa, waa_pointer_ssb, waa_pointer_totcatch, mature, percentSPR, predR, fracyr_SSB, log_SPR0, FXSPR_init);
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
  vector<Type> SPR_res_static = get_static_SPR_res(MAA, FAA_tot, which_F_age_static, pred_waa, waa_pointer_ssb, waa_pointer_totcatch, mature, percentSPR, NAA, 
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
    Type SR_a, SR_b;
    for(int y = 0; y < n_years_model + n_years_proj; y++)
    {
      log_FMSY_iter(y,0) = log(FMSY_init(y)); //starting value
      for(int a = 0; a < n_ages; a++)
      {
        M(a) = MAA(y,a);
        sel(a) = FAA_tot(y,a)/FAA_tot(y,which_F_age(y)-1); //have to look at FAA_tot to see where max F is.
        waassb(a) = pred_waa(waa_pointer_ssb-1,y,a);
        waacatch(a) = pred_waa(waa_pointer_totcatch-1, y, a);
        mat(a) = mature(y,a);
      }
      SR_a = exp(log_SR_a(y));
      SR_b = exp(log_SR_b(y));
      if(recruit_model == 3) //Beverton-Holt selected
      {
        sr_yield<Type> sryield(SR_a, SR_b, M, sel, mat, waassb, waacatch,fracyr_SSB(y),0);
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
        sr_yield<Type> sryield(SR_a, SR_b, M, sel, mat, waassb, waacatch,fracyr_SSB(y),1);
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
      log_YPR_MSY(y) = log(get_YPR(log_FMSY(y), M, sel, waacatch));
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
  REPORT(phi_mat); // NEWG
  REPORT(LAA); // NEWG
  REPORT(SDAA); // NEWG
  REPORT(GW_par); // NEWG
  REPORT(LAA_par); // NEWG
  REPORT(LW_par); // NEWG
  REPORT(pred_NAA);
  REPORT(SSB);
  REPORT(selAL);
  REPORT(selAA);
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


calibrate_cB_fn = function(
  cB_candidate = 2:9/10,
  ntrial = 1000,
  doses = c(0.05, 0.10, 0.20, 0.45, 0.65, 0.85),
  Y_B_sim = list(
    c(2.00, 2.01, 2.08, 2.76, 3.75, 4.73),
    c(2.01, 2.09, 2.24, 4.46, 5.29, 5.95),
    c(2.24, 4.00, 5.77, 5.99, 6.00, 6.00),
    c(5.04, 5.83, 5.98, 6.00, 6.00, 6.00)
  ),
  Y_T_sim = list(
    c(0.01, 0.02, 0.03, 0.06, 0.13, 0.26),
    c(0.01, 0.03, 0.04, 0.07, 0.14, 0.28),
    c(0.01, 0.02, 0.05, 0.10, 0.27, 0.55),
    c(0.01, 0.06, 0.18, 0.29, 0.51, 0.54)
  ),
  Y_R_sim = list(
    c(0.04, 0.05, 0.08, 0.20, 0.35, 0.47),
    c(0.04, 0.06, 0.09, 0.23, 0.37, 0.44),
    c(0.07, 0.14, 0.32, 0.41, 0.42, 0.44),
    c(0.22, 0.38, 0.41, 0.42, 0.44, 0.45)
  ),
  lambdaT_sim = list(
    c(0.80, 0.60, 0.60, 0.25, 0.20, 0.10),
    c(0.80, 0.40, 0.30, 0.30, 0.20, 0.35),
    c(0.40, 0.10, 0.10, 0.30, 0.30, 0.30),
    c(0.12, 0.10, 0.20, 0.30, 0.30, 0.30)
  ),
  sigma2_B_sim = c(1, 1, 1, 1),
  delta1_sim = c(3, 3, 3, 3),
  delta2_sim = c(-2, -2, -2, -2),
  delta3_sim = c(0, 0, 0, 0),
  shape_sim = c(1.5, 1.5, 1.5, 1.5),
  time_C = 24,
  target = 0.3,
  cohortsize = 3,
  ncohort = 10
){
  num_cB = length(cB_candidate)
  n_dose = length(doses)

  FP = FN = matrix(NA, ncol = num_cB, nrow = length(Y_B_sim))

  for(i_sim in 1:length(Y_B_sim)){

    Y_B = Y_B_sim[[i_sim]]
    Y_T = Y_T_sim[[i_sim]]
    Y_R = Y_R_sim[[i_sim]]
    lambdaT = lambdaT_sim[[i_sim]]

    sigma2_B = sigma2_B_sim[i_sim]
    delta1 = delta1_sim[i_sim]
    delta2 = delta2_sim[i_sim]
    delta3 = delta3_sim[i_sim]
    shape = shape_sim[i_sim]

    cutpoint = min(which(Y_B > max(Y_B) - sqrt(mean((Y_B-mean(Y_B))^2))))

    results_total = matrix(NA, ncol = n_dose, nrow = length(cB_candidate))

    for(tmp_calib in seq_along(cB_candidate)){

      monitor_cutoff_B = cB_candidate[tmp_calib]

      tau_hat_results = rep(NA, ntrial)

      for(tr in 1:ntrial){

        set.seed(tr)

        # stage I: dose-exploration
        trial = BOIN_sim(
          Y_B.true=Y_B, sigma2_B.true=sigma2_B,
          Y_T.true=Y_T,
          Y_R.true=Y_R,
          delta1.true=delta1,
          delta2.true=delta2,
          delta3.true=delta3,
          lambdaT.true=lambdaT,
          shape.true=shape,
          time_C = time_C,
          target_tox=target,
          ncohort=10, cohortsize=3, startdose=1, n.earlystop=cohortsize*ncohort,
          p.saf=0.6*target, p.tox=1.4*target, cutoff.eli=0.95, extrasafe=FALSE,
          titration=F, offset=0.05,
          max_per_dose = 9,
          monitor_cutoff_B = monitor_cutoff_B,
          seed=tr
        )

        dat_all = data.frame(d = trial$d_tol,
                             Y_B = trial$y_B_tol,
                             y = trial$y_T_tol,
                             Y_E = trial$y_R_tol,
                             Y_S = trial$y_S_tol,
                             event = trial$event_tol)

        dat_tmp = dat_all[order(dat_all$d),]
        tau_hat_results[tr] = tau_ms(dat = dat_tmp, monitor_cutoff_B = monitor_cutoff_B)

        #print(tr)
      }
      tau_hat_results = factor(tau_hat_results, levels = 1:6)
      results_total[tmp_calib, ] = as.numeric(table(tau_hat_results))
    }

    # actual = 1 but proj = 0
    FP[i_sim,] = apply(results_total[,1:(cutpoint-1),drop=FALSE], 1, sum)/ntrial
    # actual = 0 but proj = 1
    FN[i_sim,] = apply(results_total[,(cutpoint+1):n_dose,drop=FALSE], 1, sum)/ntrial
    #print(i_sim)
  }

  composite_index = sqrt((FP)^2 + (FN)^2)
  round(apply(composite_index, 2, mean), 3)

  return(list(
    composite_index = composite_index,
    best_cB = cB_candidate[which.min(round(apply(composite_index, 2, mean), 3))]
  ))
}

#calibrate_cB_fn()

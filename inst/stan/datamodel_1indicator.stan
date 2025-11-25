  // fit to the data
  // vector[N] Eta_i;
  for(i in 1:N) {
    if(held_out[i] == 0) {
       if(obs_isna[i] == 0) {
         // if (is_agg_obs[i] == 0) {
         //    eta_i[i] = inv_tr_eta(tr_eta_obs[geo_unit[i], time[i]]);
         //  } else {
         //    eta_i[i] = inv_tr_eta(tr_eta_agg[agg_unit[i], time[i]]);
         //  }
    //       y[i] ~ normal(logit(eta_i[i]), scale[i]);
    // if transformed scales for eta and dm are the same:
               y[i] ~ normal(tr_Eta_obs[geo_unit[i], time[i]], scale[i]);

      }
    }
  }


// emus
  int N_emu;
  int Ncountrytype; // no of unique ss-types
  vector[Ncountrytype] mean_log_sigma_type; // mean for each country-type
  vector<lower = 0>[Ncountrytype]  hierarchical_sigma; // sd for each country-type
  vector[N_emu] emu_roc; // observed delta emu
  vector<lower = 0>[N_emu]  sd2_emu_roc; //sd^2 of observed delta emu
  array[N_emu] int emu_for_allwomen_j; // 0 FALSE for married women pop
  array[N_emu] int c_emu_j; // geo unit index
  array[N_emu] int t_emu_j; // year index
  array[N_emu] int ss_type_index_j; // ss-type index


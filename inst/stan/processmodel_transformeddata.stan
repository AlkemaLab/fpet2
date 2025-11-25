 vector[rows(csr_extract_w(Ptilde_model_matrix))] Ptilde_model_matrix_w    = csr_extract_w(Ptilde_model_matrix);
  array[size(csr_extract_v(Ptilde_model_matrix))] int Ptilde_model_matrix_v = csr_extract_v(Ptilde_model_matrix);
  array[size(csr_extract_u(Ptilde_model_matrix))] int Ptilde_model_matrix_u = csr_extract_u(Ptilde_model_matrix);

  vector[rows(csr_extract_w(Omega_model_matrix))] Omega_model_matrix_w     = csr_extract_w(Omega_model_matrix);
  array[size(csr_extract_v(Omega_model_matrix))] int Omega_model_matrix_v  = csr_extract_v(Omega_model_matrix);
  array[size(csr_extract_u(Omega_model_matrix))] int Omega_model_matrix_u  = csr_extract_u(Omega_model_matrix);

  vector[rows(csr_extract_w(Betas_model_matrix))] Betas_model_matrix_w    = csr_extract_w(Betas_model_matrix);
  array[size(csr_extract_v(Betas_model_matrix))] int Betas_model_matrix_v = csr_extract_v(Betas_model_matrix);
  array[size(csr_extract_u(Betas_model_matrix))] int Betas_model_matrix_u = csr_extract_u(Betas_model_matrix);

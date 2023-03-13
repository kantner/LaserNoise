function [beta_sp] = func_beta_sp(par,P,N)

  % spontaneous emission (beta-)factor
  %
  % INPUT
  % - par parameter set (struct)
  % - P photon number 
  % - N carrier number   
  %
  % OUTPUT
  % - beta_sp
  
  r_sp_lasing = par.Gamma * par.vg * func_g_sp(par,P,N);
  r_sp_waste  = N.*N*par.B/par.V;
  
  beta_sp = r_sp_lasing./(r_sp_lasing + r_sp_waste);
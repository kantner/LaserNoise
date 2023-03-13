function [g_sp, dg_sp_dN] = func_g_sp_lin(par,N)

  % spontaneos emission factor without gain saturation
  %
  % INPUT
  % - P photon number 
  % - N carrier number 
  % - par parameter set (struct)
  %
  % OUTPUT
  % g_sp
  % dg_sp_dP (derivative)
  % dg_sp_dN (derivative)
  
  n    = (N/par.N_tr);
  nsq  = n.*n;
  
  g_sp     = 0.5*par.g0*log( 1 + nsq);  
  dg_sp_dN = par.g0/par.N_tr *n./( 1 + nsq);
  
end
  
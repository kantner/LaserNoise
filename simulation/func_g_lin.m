function [g, dg_dN] = func_g_lin(par,N)

  % net gain without gain saturation
  %
  % INPUT
  % - P photon number 
  % - N carrier number 
  % - par parameter set (struct)
  %
  % OUTPUT
  % g
  % dg_dP (derivative)
  % dg_dN (derivative)
  
  
  g     = par.g0*log(N./par.N_tr);  
  dg_dN = par.g0./N;
  
end
  
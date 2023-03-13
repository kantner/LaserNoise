function [g, dg_dP, dg_dN] = func_g(par,P,N)

  % net gain without gain saturation
  %
  % INPUT
  % - par parameter set (struct)
  % - P photon number 
  % - N carrier number   
  %
  % OUTPUT
  % g
  % dg_dP (derivative)
  % dg_dN (derivative)
  
  denom = (1 + par.eps*P);
  
  g     =  par.g0/denom .* log(N/par.N_tr);
  dg_dP = -par.eps*g/denom;
  dg_dN =  par.g0/denom .* (1./N);
  
  %{
  [g_lin, dg_lin_dN] = func_g_lin(par,N);
  
  g     = g_lin/denom;
  dg_dP = -par.eps*g/denom;
  dg_dN = dg_lin_dN/denom;
  %}
end
  
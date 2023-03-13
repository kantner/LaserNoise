function [g_sp, dg_sp_dP, dg_sp_dN] = func_g_sp(par,P,N)

  % spontaneos emission factor with gain saturation
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
  
  
  denom = 1 + par.eps*P;
  
  N_div_Ntr_sq = N.*N/(par.N_tr * par.N_tr);
  
  g_sp     =  0.5 * par.g0 * log(1 + N_div_Ntr_sq)/denom;
  dg_sp_dP = -par.eps*g_sp/denom;
  dg_sp_dN = par.g0 * N_div_Ntr_sq./(1 + N_div_Ntr_sq) .*1/N .*1/denom;
  
  
  %{
  denom = 1 + par.eps*P;
  
  [g_sp_lin, dg_sp_lin_dN] = func_g_sp_lin(par,N);
  
  g_sp     = g_sp_lin/denom;
  dg_sp_dP = -par.eps*g_sp/denom;
  dg_sp_dN = dg_sp_lin_dN/denom;
  %}
end
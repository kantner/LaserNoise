function [N_th, I_th] = laserThreshold(par)
  
% compute threshold carrier number and current for given parameter set

  physical_constants;
  
  N_th = par.N_tr * exp(1/(par.Gamma*par.vg*par.g0*par.tau_ph));
  I_th = elementaryCharge/par.eta * func_r(par, N_th);
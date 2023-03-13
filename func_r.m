function [r, dr_dN] = func_r(par,N)
  
  % carrier recombination rate
  %
  % INPUT
  % - par ... parameter set (struct)
  % - N   ... carrier number 
  % OUTPUT
  % - r     ... recombination rate
  % - dr_dN ... derivative

 b = par.B/par.V;
 c = par.C/(par.V*par.V);

   r    = N.*(par.A + N.*(b + N.*c));
  dr_dN = par.A + N.*(2*b + 3*c*N);
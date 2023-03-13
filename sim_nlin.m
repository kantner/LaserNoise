function [sol] = sim_nlin(par,dt,N_steps,X0,dW,F)

  % Euler-Maruyama simulation of nonlinear system
  %
  % INPUT
  %   par     ... parameters (struct)
  %   dt      ... time step size
  %   N_steps ... number of steps
  %   X0      ... initial values (P,phi,N)
  %   dW      ... white noise (ten processes)
  %   F       ... colored noise
  % OUTPUT
  %   out     ... struct
  
    physical_constants;
  
  % allocate memory
    P   = zeros(1,N_steps);
    phi = zeros(1,N_steps);
    N   = zeros(1,N_steps);
  
  % initial values
    P(1)   = X0(1);
    phi(1) = X0(2);
    N(1)   = X0(3);
  
  % aux
    Gamma_vg  = par.Gamma * par.vg;
    injection = par.eta*par.I/elementaryCharge;
    sqrt_inj  = sqrt(injection);
    
  % CW frequency  
    Omega = zeros(1,N_steps);
    
  for i = 1 : N_steps-1
    % gain
      Gamma_vg_g     = Gamma_vg * func_g(par,P(i),N(i));
      Gamma_vg_g_sp  = Gamma_vg * func_g_sp(par,P(i),N(i));
      Gamma_vg_g_abs = Gamma_vg_g_sp - Gamma_vg_g;
    
    % reco
      r = func_r(par,N(i));      
      
    % noise due to spontaneous and stimulated recombination
      light_matter_noise = + sqrt(Gamma_vg_g_sp * P(i) ) * dW(5,i) ...
                           + sqrt(Gamma_vg_g_abs * P(i)) * dW(7,i) ...
                           + sqrt(Gamma_vg_g_sp) * dW(9,i);
      
      
    % update  
      P(i+1)   = P(i) + ...
                 + ( -(P(i) - par.P_th)/par.tau_ph  +  Gamma_vg_g * P(i) + Gamma_vg_g_sp + par.sigma(1)*2*P(i) * F(1,i) ) * dt ...
                 + sqrt(P(i)*(1+par.P_th)/par.tau_ph) * dW(1,i) ...
                 + sqrt(par.P_th*(1+P(i))/par.tau_ph) * dW(3,i) ...
                 + light_matter_noise;
                 
      Omega(i+1) = 0.5*par.alpha_H * Gamma_vg_g;
      
      phi(i+1) = phi(i) + (Omega(i+1) + par.sigma(2)*F(2,i)) * dt ...
                 + 0.5/P(i) * (  sqrt(P(i)*(1+par.P_th)/par.tau_ph) * dW(2,i) ...
                               + sqrt(par.P_th*(1+P(i))/par.tau_ph) * dW(4,i) ...
                               + sqrt(Gamma_vg_g_sp  * P(i)) * dW(6,i) ...
                               + sqrt(Gamma_vg_g_abs * P(i)) * dW(8,i) ...
                               + sqrt(Gamma_vg_g_sp) * dW(10,i) ...
                               );

      N(i+1)   = N(i) + (injection - r - Gamma_vg_g * P(i) - Gamma_vg_g_sp + par.sigma(3) * sqrt(N(i)) * F(3,i) ) * dt ...
                 + sqrt_inj * dW(11,i) ...
                 + sqrt(r) * dW(12,i) ...
                 - light_matter_noise;
      
  end
  
  sol.P     = P;
  sol.phi   = phi;
  sol.N     = N;
  sol.Omega = Omega;
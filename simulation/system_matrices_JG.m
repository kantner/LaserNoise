function [J,G] = system_matrices_JG(par)
 
    physical_constants;
  
  % matrices for system [P, phi, N]  
    
  % steady state
    [P_ss, N_ss] = steady_state(par);

  % aux
    Gamma_vg  = par.Gamma * par.vg;
    injection = par.eta*par.I/elementaryCharge;
    
    [r, dr_dN] = func_r(par,N_ss); 
    
    [g,    dg_dP,    dg_dN]    = func_g(par,P_ss,N_ss);
    [g_sp, dg_sp_dP, dg_sp_dN] = func_g_sp(par,P_ss,N_ss);
     g_abs                     = g_sp - g;
    
  % relaxation rates   
    Gamma_N = dr_dN    + Gamma_vg*(dg_dN.*P_ss + dg_sp_dN);
    Gamma_P = 1/par.tau_ph - Gamma_vg*(g + P_ss.*dg_dP + dg_sp_dP);
     
    
  % Jacobian
    J = [ -Gamma_P,                             0,  Gamma_N - dr_dN;
          0.5*par.alpha_H*Gamma_vg*dg_dP,       0,  0.5*par.alpha_H*Gamma_vg*dg_dN;
          -Gamma_vg*(g+dg_dP*P_ss + dg_sp_dP),  0,  -Gamma_N];
  
  % noise matrix    
    G = [sqrt((1+par.P_th)*P_ss/par.tau_ph), 0,                                             sqrt(par.P_th*(1+P_ss)/par.tau_ph), 0,                                              sqrt(Gamma_vg*g_sp*P_ss),   0,                              sqrt(Gamma_vg*g_abs*P_ss),  0,                              sqrt(Gamma_vg*g_sp),  0,                                    0,               0;
         0,                                  sqrt((1+par.P_th)*P_ss/(4*P_ss^2*par.tau_ph)), 0,                                  sqrt(par.P_th*(1+P_ss)/(4*P_ss^2*par.tau_ph)),  0,                          sqrt(0.25*Gamma_vg*g_sp/P_ss),  0,                          sqrt(0.25*Gamma_vg*g_abs/P_ss), 0,                    sqrt(0.25*Gamma_vg*g_sp/(P_ss*P_ss)), 0,               0;
         0,                                  0,                                             0,                                  0,                                             -sqrt(Gamma_vg*g_sp*P_ss),   0,                             -sqrt(Gamma_vg*g_abs*P_ss),  0,                             -sqrt(Gamma_vg*g_sp),  0,                                    sqrt(injection), sqrt(r)];
     
     
    J = sparse(J); 
    G = sparse(G);
    
    

  return
  % relaxation rates
    

     
     Gamma_R = 0.5 * (Gamma_P + Gamma_N)
     Omega_R = sqrt( Gamma_vg*(dg_dN .* P_ss + dg_sp_dN).*Gamma_vg.*(g + dg_dP.* P_ss + dg_sp_dP)  - 0.25*(Gamma_N - Gamma_P)^2 )
    
     eigs(J)
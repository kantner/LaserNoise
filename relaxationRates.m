function [Gamma_R, Omega_R, Gamma_P, Gamma_N] = relaxationRates(par)


  % matrices for system [P, phi, N]  
    
  % steady state
    [P_ss, N_ss] = steady_state(par);

  % aux
    Gamma_vg  = par.Gamma * par.vg;
    
    [r, dr_dN] = func_r(par,N_ss); 
    
    [g,    dg_dP,    dg_dN]    = func_g(par,P_ss,N_ss);
    [g_sp, dg_sp_dP, dg_sp_dN] = func_g_sp(par,P_ss,N_ss);

  % rates  
    Gamma_P = 1/par.tau_ph - Gamma_vg*(g + P_ss.*dg_dP + dg_sp_dP);
    Gamma_N = dr_dN    + Gamma_vg*(dg_dN.*P_ss + dg_sp_dN);
        
    %Gamma_R = 0.5 * (Gamma_P + Gamma_N)
    %Omega_R = real(sqrt( Gamma_vg*(dg_dN .* P_ss + dg_sp_dN).*Gamma_vg.*(g + dg_dP.* P_ss + dg_sp_dP)  - 0.25*(Gamma_N - Gamma_P)^2 ));
    
  % Jacobian
    [J,G] = system_matrices_JG(par);
    
    lambda = eigs(J);
    
  % delete trivial eigenvalue
    idx = find(lambda==0);
    lambda(idx)=[];
    
    
    if abs(imag(lambda)) > 0
       % complex conjugate pair
         Gamma_R = -real(lambda(1));
         Omega_R = abs(imag(lambda(1)));
    else
       % two real eigenvalues
         Gamma_R = 0; % we do not assign a relaxation rate in this case
         Omega_R = 0;
       
       
    end
    
    
     
function [P_ss, N_ss] = steady_state(par)

  
  % compute steady state of noise-free system by solving a scalar condition
    
  % load physical constants   
    physical_constants;
    
     
     
  %% debug mode (0 = off, 1 = on)
     DEBUG = 0;
  
    
  %% P as function of N  
      func_P = @(N)  par.P_th + par.tau_ph * ( par.eta*par.I/elementaryCharge - func_r(par, N) );
     %dfunc_P = @(N) -par.tau_ph * ( par.A  + N.*(2*par.B/par.V  + 3*par.C/(par.V^2)* N));
   

      
     
     
  %% scalar equation
     f = @(N) func_g_lin(par, N)/par.g0 .* func_P(N) ...
              + func_g_sp_lin(par, N)/par.g0 ...
              - (1 + par.eps*func_P(N)).*(func_P(N)-par.P_th)./(par.Gamma*par.vg*par.g0*par.tau_ph);
          
     %df = @(N) (func_P(N) + 1./(1 + (par.N_tr./N).* (par.N_tr./N)))./N + (log(N./par.N_tr) - (1 + 2*par.eps*func_P(N) )/(par.Gamma*par.vg*par.g0*par.tau_ph)).*dfunc_P(N);
     
     
     
         
  %% determine lower and upper bounds for N for bisection

  if par.eps > 0 % case with nonlinear gain
      
     % lower bound N_min [ using g <= g_sp <= 0.5*x^2 and (1+P)/(1+eps P) <= max(1,1/eps) ] - 3rd order polynomial
        a0 = -par.eta*par.I/(elementaryCharge*par.Gamma*par.vg*par.g0);
        a1 = par.A*par.N_tr/(par.Gamma*par.vg*par.g0);
        a2 = 0.5 * max([1,1/par.eps]) +  par.B*(par.N_tr^2)/(par.V*par.Gamma*par.vg*par.g0);
        a3 = par.C*(par.N_tr^3)/((par.V^2)*par.Gamma*par.vg*par.g0);
     
        %X_low    = zeros(3,1);
        %X_low(1) = -(a2/(3*a3)) - (2^(1/3)*(-a2^2 + 3*a1*a3))/(3*a3*(-2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2 + sqrt(4*(-a2^2 + 3*a1*a3)^3 + (-2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2)^2))^(1/3)) + (-2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2 + sqrt(4*(-a2^2 + 3*a1*a3)^3 + (-2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2)^2))^(1/3)/ (3*2^(1/3)*a3);
        %X_low(2) = -(a2/(3*a3)) + ((1 + 1i*sqrt(3))*(-a2^2 + 3*a1*a3))/ (3*2^(2/3)*a3*(-2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2 + sqrt(4*(-a2^2 + 3*a1*a3)^3 + (-2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2)^ 2))^(1/3)) - ((1 - 1i*sqrt(3))*(-2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2 + sqrt(4*(-a2^2 + 3*a1*a3)^3 + (-2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2)^2))^(1/3))/(6*2^(1/3)*a3);
        %X_low(3) = -(a2/(3*a3)) + ((1 - 1i*sqrt(3))*(-a2^2 + 3*a1*a3))/(3*2^(2/3)*a3* (-2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2 + sqrt(4*(-a2^2 + 3*a1*a3)^3 + (-2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2)^2))^(1/3)) - ((1 + 1i*sqrt(3))*(-2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2 + sqrt(4*(-a2^2 + 3*a1*a3)^3 + (-2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2)^2))^(1/3))/(6*2^(1/3)*a3);
        %X_low
        X_low = roots([a3 a2 a1 a0]);
        %f_low(X_low*par.N_tr)
        % select
        N_min = par.N_tr * real(X_low(real(X_low)>0));
        
      % upper bound N_max [ using (x-1)/x <= g <= g_sp and (1+P)/(1+eps P) >= min(1,1/eps) ]
        b0 = -min([1,1/par.eps]);
        b1 = min([1,1/par.eps]) - par.eta*par.I/(elementaryCharge*par.Gamma*par.vg*par.g0);       
        b2 = par.A*par.N_tr/(par.Gamma*par.vg*par.g0);
        b3 = par.B*(par.N_tr^2)/(par.V*par.Gamma*par.vg*par.g0);
        b4 = par.C*(par.N_tr^3)/((par.V)^2*par.Gamma*par.vg*par.g0);
  
        %X_up    = zeros(4,1);
        %X_up(1) = -(b3/(4*b4)) - (1/2)*sqrt(b3^2/(4*b4^2) - (2*b2)/(3*b4) + (2^(1/3)*(b2^2 - 3*b1*b3 + 12*b0*b4))/ (3*b4*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^(1/3)) + (1/(3*2^(1/3)*b4))*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3)) - (1/2)*sqrt(b3^2/(2*b4^2) - (4*b2)/(3*b4) - (2^(1/3)*(b2^2 - 3*b1*b3 + 12*b0*b4))/ (3*b4*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^(1/3)) - (1/(3*2^(1/3)*b4))*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3) - (-(b3^3/b4^3) + (4*b2*b3)/b4^2 - (8*b1)/b4)/(4*sqrt(b3^2/(4*b4^2) - (2*b2)/(3*b4) + (2^(1/3)*(b2^2 - 3*b1*b3 + 12*b0*b4))/(3*b4*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3)) + (1/(3*2^(1/3)*b4))*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt( -4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3))));
        %X_up(2) = -(b3/(4*b4)) - (1/2)*sqrt(b3^2/(4*b4^2) - (2*b2)/(3*b4) + (2^(1/3)*(b2^2 - 3*b1*b3 + 12*b0*b4))/ (3*b4*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^(1/3)) + (1/(3*2^(1/3)*b4))*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3)) + (1/2)*sqrt(b3^2/(2*b4^2) - (4*b2)/(3*b4) - (2^(1/3)*(b2^2 - 3*b1*b3 + 12*b0*b4))/ (3*b4*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^(1/3)) - (1/(3*2^(1/3)*b4))*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3) - (-(b3^3/b4^3) + (4*b2*b3)/b4^2 - (8*b1)/b4)/(4*sqrt(b3^2/(4*b4^2) - (2*b2)/(3*b4) + (2^(1/3)*(b2^2 - 3*b1*b3 + 12*b0*b4))/(3*b4*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3)) + (1/(3*2^(1/3)*b4))*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt( -4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3))));
        %X_up(3) = -(b3/(4*b4)) + (1/2)*sqrt(b3^2/(4*b4^2) - (2*b2)/(3*b4) + (2^(1/3)*(b2^2 - 3*b1*b3 + 12*b0*b4))/ (3*b4*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^(1/3)) + (1/(3*2^(1/3)*b4))*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3)) - (1/2)*sqrt(b3^2/(2*b4^2) - (4*b2)/(3*b4) - (2^(1/3)*(b2^2 - 3*b1*b3 + 12*b0*b4))/ (3*b4*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^(1/3)) - (1/(3*2^(1/3)*b4))*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3) + (-(b3^3/b4^3) + (4*b2*b3)/b4^2 - (8*b1)/b4)/(4*sqrt(b3^2/(4*b4^2) - (2*b2)/(3*b4) + (2^(1/3)*(b2^2 - 3*b1*b3 + 12*b0*b4))/(3*b4*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3)) + (1/(3*2^(1/3)*b4))*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt( -4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3))));
        %X_up(4) = -(b3/(4*b4)) + (1/2)*sqrt(b3^2/(4*b4^2) - (2*b2)/(3*b4) + (2^(1/3)*(b2^2 - 3*b1*b3 + 12*b0*b4))/ (3*b4*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^(1/3)) + (1/(3*2^(1/3)*b4))*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3)) + (1/2)*sqrt(b3^2/(2*b4^2) - (4*b2)/(3*b4) - (2^(1/3)*(b2^2 - 3*b1*b3 + 12*b0*b4))/ (3*b4*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^(1/3)) - (1/(3*2^(1/3)*b4))*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3) + (-(b3^3/b4^3) + (4*b2*b3)/b4^2 - (8*b1)/b4)/(4*sqrt(b3^2/(4*b4^2) - (2*b2)/(3*b4) + (2^(1/3)*(b2^2 - 3*b1*b3 + 12*b0*b4))/(3*b4*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt(-4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3)) + (1/(3*2^(1/3)*b4))*(2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4 + sqrt( -4*(b2^2 - 3*b1*b3 + 12*b0*b4)^3 + (2*b2^3 - 9*b1*b2*b3 + 27*b0*b3^2 + 27*b1^2*b4 - 72*b0*b2*b4)^2))^ (1/3))));
        X_up = roots([b4 b3 b2 b1 b0]);
        % select
        N_max = par.N_tr * real(X_up(real(X_up)>0));
        %X_up
        %f_up(X_up*par.N_tr)    
        
        
        if DEBUG == 1
           
           a0
           a1
           a2
           a3
           X_low
           N_min
           roots([a3 a2 a1 a0])
           
           
           
           
           b0
           b1
           b2
           b3
           b4           
           X_up
           N_max
           roots([b4 b3 b2 b1 b0])
            
           N_range = logspace(log10(N_min),log10(N_max),10001);
           figure(99999);clf;hold all;
           plot(N_range,f(N_range),'r-')
           plot(N_range,0*N_range,'k--')
            
            
        end
        
        
  elseif par.eps == 0
     
    % lower bound: 5th order polynomial  
      a5 = -0.5*par.tau_ph * par.C *par.N_tr^3/(par.V^2);
      a4 = -0.5*par.tau_ph * par.B *par.N_tr^2/par.V;
      a3 = -0.5*par.tau_ph * par.A *par.N_tr +  1/(par.Gamma*par.vg*par.g0*par.tau_ph) *  (par.tau_ph * par.C *par.N_tr^3/(par.V^2));
      a2 = 0.5*par.P_th + 0.5*par.tau_ph*par.eta*par.I/elementaryCharge + 0.5 + 1/(par.Gamma*par.vg*par.g0*par.tau_ph) *  (par.tau_ph * par.B *par.N_tr^2/par.V);
      a1 = 1/(par.Gamma*par.vg*par.g0*par.tau_ph) *  (par.tau_ph * par.A *par.N_tr);
      a0 = -1/(par.Gamma*par.vg*par.g0*par.tau_ph) *  (par.tau_ph * par.eta * par.I/elementaryCharge);
      
      X_low= roots([a5 a4 a3 a2 a1 a0]);
      N_min = par.N_tr * real(X_low(real(X_low)>0));
      N_min = min(N_min);

      
    % upper bound: threshold density
      %b4 = (1/(par.Gamma*par.vg*par.g0*par.tau_ph) - 1) * par.tau_ph * par.C *par.N_tr^3/(par.V^2);
      %b3 = (1/(par.Gamma*par.vg*par.g0*par.tau_ph) - 1) * par.tau_ph * par.B *par.N_tr^2/par.V   + par.tau_ph * par.C *par.N_tr^3/(par.V^2);
      %b2 = (1/(par.Gamma*par.vg*par.g0*par.tau_ph) - 1) * par.tau_ph * par.A *par.N_tr   + par.tau_ph * par.B *par.N_tr^2/par.V;
      %b1 = (1/(par.Gamma*par.vg*par.g0*par.tau_ph) - 1) * par.tau_ph *par.eta*par.I/elementaryCharge + par.P_th + 1  + par.tau_ph * par.A *par.N_tr;
      %b0 = -par.P_th - 1  - par.tau_ph * par.eta*par.I/elementaryCharge;
            
      %X_up  = roots([b4 b3 b2 b1 b0])
      %N_max = par.N_tr * real(X_up(real(X_up)>0))
      %N_max = min(N_max);

      N_th = laserThreshold(par);
      N_max = N_th;
      
      if DEBUG == 1
           
          
           N_range = logspace(log10(N_min),log10(N_max),10001);
           N_range = logspace(-5,10,100001);
           figure(99999);clf;hold all;
           plot(N_range,f(N_range),'r-')
           plot(N_range,0*N_range,'k--')
           
           plot([1 1]*N_min,[-1 1]*1E6,'b--')
           plot([1 1]*N_max,[-1 1]*1E6,'r--')
           plot([1 1]*par.N_tr,[-1 1]*1E6,'k--')
           set(gca,'XScale','log')
           ylim([-1 1]*1E6)
            
            
      end
        
      
      
  end
        
        
   %% solve        
   
      if sign(f(N_min)) == sign(f(N_max))
          
        fprintf(1,'steadyState.m: Interval endpoints do not differ in sign ... return NaN\n')
        N_ss = NaN;
        P_ss = NaN;
          
      else
   
   
        opts = optimset('Display','none','TolX',1E-36);        
        [N_ss,~,EXITFLAG] = fzero(f,[N_min, N_max],opts);      
        P_ss = func_P(N_ss);
        
        if DEBUG == 1
            EXITFLAG
            N_ss
            P_ss
            
        end 
        
        
      % EXITFLAG
      %  1  fzero found a zero X.
      % -1  Algorithm terminated by output function.
      % -3  NaN or Inf function value encountered during search for an interval containing a sign change.  
      % -4  Complex function value encountered during search for an interval containing a sign change.
      % -5  fzero may have converged to a singular point.
      % -6  fzero can not detect a change in sign of the function.
        if EXITFLAG ~= 1
          N_ss = NaN;
          P_ss = NaN;
        end
      
      end
      
      
      

    
      
      
      
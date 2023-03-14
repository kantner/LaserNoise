function [out] = sim_lin(par,dt,N_steps,X0,dW,F)

  % Euler-Maruyama simulation of linearized system
  %
  % INPUT
  %   par     ... parameters (struct)
  %   dt      ... time step size
  %   N_steps ... number of steps
  %   X0      ... initial values (P,phi,N)
  %   dW      ... white noise
  %   F       ... colored noise
  % OUTPUT
  %   out     ... struct
  
  
  % allocate memory
    X   = zeros(3,N_steps);
    
  % initial values
    X(1,1) = X0(1);
    X(2,1) = X0(2);
    X(3,1) = X0(3);
     
  % generate matrices  
    [J, G] = system_matrices_JG(par);
  
  % steady state values  
    [P_ss, N_ss] = steady_state(par);
    
    A = speye(3) + J*dt;
    
  % iterate  
    for i = 1 : N_steps-1
      X(:,i+1) =    A * X(:,i) ...                                       % deterministic
                  + G * dW(:,i) ...                                      % white noise
                  + [1;0;0] * par.sigma(1) * 2*P_ss     * F(1,i) * dt ...  % colored noise P
                  + [0;1;0] * par.sigma(2)              * F(2,i) * dt ...  % colored noise phi
                  + [0;0;1] * par.sigma(3) * sqrt(N_ss) * F(3,i) * dt;     % colored noise N
                
    end
  
  % output  
    out.P   = X(1,:);
    out.phi = X(2,:);
    out.N   = X(3,:);
  
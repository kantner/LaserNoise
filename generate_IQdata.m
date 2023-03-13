function [I, Q, I_no_noise, Q_no_noise] = generate_IQdata(P, phi, l, dt, deltaOmega, eta_det, eta_out, sigma_IQ)

% generate I and Q time series of self-delayed heterodyne beat note measurement
% from simulated time series
% INPUT
%   P   .......... power time series    (P(t)   = P_ss       + dP(t)   )
%   phi .......... phase time series    (phi(t) = Omega_CW*t + dphi(t) )
%   l   .......... discrete delay time
%   dt  .......... time step
%   deltaOmega ... frequency offset (AOM etc.)
%   eta_det ...... detector efficiencity
%   sigma ........ measurement noise (standard deviation)
% OUTPUT
%   I ............ time series inlcuding simulated measurement noise
%   Q ............ time series inlcuding simulated measurement noise
%   I_no_noise ... time series without simulated measurement noise
%   Q_no_noise.... time series without simulated measurement noise

N  = length(P);
t  = [0:N-1]*dt;

ampl = eta_det * eta_out * sqrt( P .* circshift(P, l));
arg  = phi - circshift(phi, l) - deltaOmega * t;

I_no_noise =  ampl .* cos(arg);
Q_no_noise =  ampl .* sin(arg);


I_noise = sigma_IQ * sqrt(dt) * randn(1,N);
Q_noise = sigma_IQ * sqrt(dt) * randn(1,N); % uncorrelated (!?)

I = I_no_noise + I_noise;
Q = Q_no_noise + Q_noise;




% delete leading l values   
  I(1:l)          = [];
  Q(1:l)          = [];
  
  I_no_noise(1:l) = [];
  Q_no_noise(1:l) = [];
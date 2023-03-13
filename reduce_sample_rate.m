function [x_coarse, t_coarse] = reduce_sample_rate(x, t, n, mode)
% Reduce sampling rate of time series by coarse-graining.
%
% INPUT
%   x ...... time series
%   t ...... column vector of sampling times
%   n ...... integer that determines coarse graining (n=1 gives back original time series)
%   mode ... frequency space truncation (mode = 1), averaging (mode = 2) or stroboscopic (mode = 3)
% OUTPUT
%   x_coarse ...... coarse grained time series
%   t_coarse ...... vector of coarse grained sampling times
%
% The time series is required to have a length of N = 2^Integer, otherwise it will be
% trunctated. Similarly, the coarse graining index n = 2^Integer is required to be a
% power of 2.


% check if x is a column vector
  if size(x,1) == 1 && size(x,2) > 1
    % is column vector
  else
    error('x is required to be a column vector.')
  end
  
% check if n is a power of 2
  if log2(n) ~= round(log2(n))
    error('n is required to be a power of 2: n = 2^Integer.') 
  end
  
% truncate x
  N = 2^floor(log2(length(x)));
  x = x(1:N);
  t = t(1:N);
  dt = t(2)-t(1);
  
  f = [-N/2 : N/2-1]/(N*dt);
  
% coarse graining
  dt_coarse = dt * n;
  N_coarse  = N/n;
  
  t_coarse  = [0 : N_coarse-1]*dt_coarse;
  
  
  switch(mode)
    case 1  % frequency space truncation
      f_coarse  = [-N_coarse/2 : N_coarse/2-1]/(N_coarse*dt_coarse);  
      
      idx1 = find(f_coarse(1)   == f);
      idx2 = find(f_coarse(end) == f);
      
      X = fftshift(fft(x));
      
      X = X(idx1 : idx2);
      x_coarse = 1/n * real(ifft(ifftshift(X)));
              
    case 2 % averaging
      x_coarse = mean(reshape(x,n,N_coarse),1);        
    case 3 % stroboscopic
      x_coarse = x(1:n:end);
  end
  
  
  
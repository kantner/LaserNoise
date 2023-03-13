function [w] = hann_window(varargin)
% Replica of the MATLAB function hann.m (part of the Signal Processing Toolbox)
% 
% - hann_window(N) returns the N-point symmetric Hann window in a column vector. 
% - hann_window(N,SFLAG) generates the N-point Hann window using SFLAG window sampling.
%   SFLAG may be either 'symmetric' or 'periodic'. By default, a symmetric
%   window is returned. 



  narginchk(1, 2)

  N = varargin{1};
  if nargin == 1
    % set sflag to 'symmetric' (MATLAB default)
      sflag = 'symmetric';
  elseif nargin == 2
      sflag = varargin{2};
  end
  
  % create symmetric Hann window with N+1 points  
    w = 0.5*(1 - cos(2*pi*[0:N]'/N));
  
  if strcmp(sflag, 'periodic')
     % delete last point
       w(end) = [];
  elseif strcmp(sflag, 'symmetric')
     % leave it as it is
  else
     % error
       sflag
       error('undefined option')       
  end
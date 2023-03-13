function [S] = PSD(varargin)
% estimate power spectral density


  narginchk(1, 3) % min 1 input, max 3 inputs

  x  = varargin{1};
  if nargin == 1
    % input: x ... set others to default
      dt                 = 1;
      hann_window_switch = 0;
  elseif nargin == 2
    % input: x, dt ... set others to default
      dt                 = varargin{2};
      hann_window_switch = 0;
  elseif nargin == 3
    % input: x, dt, hann_window_switch
      dt                 = varargin{2};
      hann_window_switch = varargin{3};
  end
  
  switch(hann_window_switch)
    case 0 % Hann window off
      S = dt/length(x) * abs(fftshift(fft(x ))).^2;
    case 1 % Hann window on
      S = dt/length(x) * abs(fftshift(fft(x .* hann_window(length(x),'periodic')' ))).^2;
      %PSD = @(x, dt) dt/length(x) * abs(fftshift(fft(x .* hann(length(x),'periodic')' ))).^2;  <<------ this version requires MATLAB's Signal Processing Toolbox (separate license)
    otherwise
      error('not defined')  
  end
   
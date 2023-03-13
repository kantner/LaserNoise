function [W]=fbm1d(H,n)
% fast one dimensional fractional Brownian motion (FBM) generator
% output is 'W_t' with t in [0,T] using 'n' equally spaced grid points;
% code uses Fast Fourier Transform (FFT) for speed.
% INPUT:
%          - Hurst parameter 'H' in [0,1]
%          - number of grid points 'n', where 'n' is a power of 2;
%            if the 'n' supplied is not a power of two,
%            then we set n=2^ceil(log2(n)); default is n=2^12;
%          - final time 'T'; default value is T=1;
% OUTPUT:
%          - Fractional Brownian motion 'W_t' for 't';
%          - time 't' at which FBM is computed;
%            If no output it invoked, then function plots the FBM.
% Example: plot FBM with hurst parameter 0.95 on the interval [0,10]
% [W,t]=fbm1d(0.95,2^12,10); plot(t,W)

% Reference: 
% Kroese, D. P., & Botev, Z. I. (2015). Spatial Process Simulation. 
% In Stochastic Geometry, Spatial Statistics and Random Fields(pp. 369-404)
% Springer International Publishing, DOI: 10.1007/978-3-319-10064-7_12

%if (H>1)|(H<0) % Hurst parameter error check
%    error('Hurst parameter must be between 0 and 1')
%end
if nargin<2
    n=2^12;  % grid points
else
    n=2^ceil(log2(n));
end

VH = 1/(sin(pi*H)*gamma(2*H+1));

r=nan(n+1,1);
r(1) = 1;
idx=1:n;
r(idx+1) = VH * 0.5*((idx+1).^(2*H) - 2*idx.^(2*H) + (idx-1).^(2*H));
r=[r; r(end-1:-1:2)]; % first row of circulant matrix
lambda=real(fft(r))/(2*n); % eigenvalues
W=fft(sqrt(lambda).*complex(randn(2*n,1),randn(2*n,1)));
W = n^(-H)*cumsum(real(W(1:n))); % rescale


if nargout==0
    plot(W); title('Fractional Brownian motion');
    xlabel('time $t$','interpreter','latex')
    ylabel('$W_t$','interpreter','latex')
end

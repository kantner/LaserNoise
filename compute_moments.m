function [mu, Sigma] = fit_gaussian(logL, x)
%
% Fit to multi-variate Gaussian
% p(x, mu, Sigma) = 1/sqrt(det(2 pi Sigma)) * exp(-0.5 (x-mu)' inv(Simga) (x-mu) )
%
% Method: Computation of moments
%
% INPUT
%   logL ... log-likelihood data (N x 1)
%   x ...... parameter vector (N x m)

% OUTPUT
%   mu
%   Sigma



  m = size(x,2);

  
  norm = sum(logL);
  
  
% mean
  mu = zeros(1,m);
  for i = 1 : m
      mu(i) = sum(logL.*x(:,i))/norm;
  end
  
  
% cov
  Sigma = zeros(m,m);
  for i = 1 : m
    for j = i : m
      Sigma(i,j) = sum(logL.*(x(:,i) - mu(i)).*(x(:,j)-mu(j)) )/norm;
      if j > i
      Sigma(j,i) = Sigma(i,j);
      end
    end
  end
  
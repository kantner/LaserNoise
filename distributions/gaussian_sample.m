function [X] = gaussian_sample(mu,sigmasq,n)
  
  X = mu + sqrt(sigmasq)*randn(n,1);

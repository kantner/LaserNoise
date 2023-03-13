function [X] = lognormal_sample(mu,sigmasq,n)

  X = exp(mu + sqrt(2*sigmasq) * erfinv(2*rand(n,1)-1));
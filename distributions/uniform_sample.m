function [X] = uniform_sample(xmin,xmax,n)
  
  X = xmin + rand(n,1)*(xmax-xmin);

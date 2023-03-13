function [p] = uniform_pdf(x,xmin,xmax)
  p = heaviside(x-xmin).*heaviside(xmax-x)./(xmax-xmin);

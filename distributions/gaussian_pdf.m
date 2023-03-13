function [p] = gaussian_pdf(x,mu,sigmasq)

p = exp(-0.5*(x-mu).*(x-mu)/sigmasq)/sqrt(2*pi*sigmasq);

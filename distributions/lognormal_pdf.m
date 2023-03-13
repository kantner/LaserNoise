function [p] = lognormal_pdf(x,mu,sigmasq)

  tmp = log(x)-mu;

  p = exp(-0.5*tmp.*tmp/sigmasq)./(x*sqrt(2*pi*sigmasq));

  p(x<=0) = 0;
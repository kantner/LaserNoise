function [logL] = loglikelihood_Sxx(x, f, S_xx)
  
    % compute log-likelihood for S_zz PSD
    %
    % INPUT
    %   x ............. parameter vector x = [C nu S_inf]
    %   f ............. vector of frequencies
    %   S_xx .......... data at respective frequencies
    %
    % OUTPUT
    %   logL .......... log-likelihood for given data and parameters
    
    
    
    
   %% model PSDs for inference 
    % signal
      S_xx_model   = x(1)./(f.^x(2)) + x(3);
      
   %% statistical model for S_xx data 
    % Assume exponentially distributed data (chi^2_{k=2})
    %
    %   p(x) =  lambda .* exp(-lambda.*x)
    %
    % where x is the data and the parameter lambda controls the moments of
    % the distribution. For the mean value it holds
    %
    %   E(x) = int_0^\infty dx x p(x) = 1/lambda.
    %
    % Here the lambda is frequency-dependent (different at each data point)
    % and always identified with the inverse of the modeled (!) data:
    %  
    %   lambda =  1/S_xx
    %
     
    % set vector of lambdas
      lambda = 1./S_xx_model;
      
    % Instead of p, we comoute log(p) directly, where x is identified with observed (!) data
      logp = log(lambda) - lambda.*S_xx;
        
    % compute log-likelihood  
      logL = sum(logp);
           
    
end
 
  
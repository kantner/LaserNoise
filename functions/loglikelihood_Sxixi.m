function [logL] = loglikelihood_Sxixi(x, f, S_xixi)     

    % compute log-likelihood for S_xixi PSD
    %
    % INPUT
    %   x ............. parameter \sigma
    %   f ............. vector of frequencies
    %   S_xixi ........ data at respective frequencies
    %
    % OUTPUT
    %   logL .......... log-likelihood for given data and parameters
    
    
    
    
   %% model PSD for inference 
      S_xixi_model = x^2 * f.^2;
    
      
   %% statistical model for S_xixi data 
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
    %   lambda =  1/S_xixi
    %
     
    % set vector of lambdas
      lambda = 1./S_xixi_model;
      
    % Instead of p, we comoute log(p) directly, where x is identified with observed (!) data
      logp = log(lambda) - lambda.*S_xixi;
        
    % compute log-likelihood  
      logL = sum(logp);
      
end
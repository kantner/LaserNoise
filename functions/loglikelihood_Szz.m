function [logL] = loglikelihood_Szz(x, f, S_zz, Hsq, plot_switch)
  
    % compute log-likelihood for S_zz PSD
    %
    % INPUT
    %   x ............. parameter vector x = [C nu S_inf sigma]
    %   f ............. vector of frequencies
    %   S_zz .......... data at respective frequencies
    %   Hsq ........... absolute square of transfer function at respective frequencies
    %   plot_switch ... show plot of data S_zz and fits with chosen parameters x
    %
    % OUTPUT
    %   logL .......... log-likelihood for given data and parameters
    
    
    
    
   %% model PSDs for inference 
    % signal
      S_xx_model   = x(1)./(f.^x(2)) + x(3);
      
    % noise  
    
      S_xixi_model = x(4)^2 * f.^2;
    % raw data  
      S_zz_model   = Hsq.*S_xx_model + S_xixi_model;
      
   %% statistical model for S_zz data 
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
    %   lambda =  1/S_zz
    %
     
    % set vector of lambdas
      lambda = 1./S_zz_model;
      
    % Instead of p, we comoute log(p) directly, where x is identified with observed (!) data
      logp = log(lambda) - lambda.*S_zz;
        
    % compute log-likelihood  
      logL = sum(logp);
      
      
      if plot_switch == 1
        figure(1000);clf;hold all;
          plot(f, S_zz,'b-','LineWidth',1,'DisplayName','S_{z,z} (data)')
          plot(f, S_xx_model,'k-','LineWidth',3,'DisplayName','S_{x,x} (model)')
          plot(f, S_xixi_model,'m-','LineWidth',2,'DisplayName','S_{\xi,\xi} (model)')
          plot(f, S_zz_model,'c-','LineWidth',2,'DisplayName','S_{z,z} (model)')
          set(gca,'XScale','log')
          set(gca,'YScale','log')
          box on
          title(['parameters: C = ',num2str(x(1),'%g'),...
                 ' | \nu = ',num2str(x(2),'%g'),...
                 ' | S_{\infty} = ',num2str(x(3),'%g'),...
                 ' | \sigma = ',num2str(x(4),'%g')])                 
          xlabel('frequency (Hz)')
          ylabel('S_{z,z}')
          axis tight
          drawnow
      end
    
  end
  
  
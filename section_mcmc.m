  fprintf(1,'\n') 
  fprintf(1,'=== filtering of SDH signal ===\n')
  % The reconstruction is for the FN-PSD only.
  
  %% convolution function   
   % convolution function h(t) and Fourier transform H(\omega) for phase signal       
     h                   = zeros(N_meas,1);
     h(1)                = 1;  % instantaneous signal
     h(N_meas+1-l_meas)  = -1; % retarded signal
     
     H                   = transpose(fftshift(fft(h)));
     
     Hsq = abs(H).^2;
     
        
  %% hidden signal and noise
     S_xx   = (2*pi*f).^2 .* PSD(sol_nlin.phi - Omega_CW*t, dt, par.hann_window_switch);
     S_xixi = (2*pi*f_meas).^2 .* PSD(phase_noise, dt_meas, par.hann_window_switch);     
  
   % reference value from hidden data
     idx = 1E7 <= f & f<=5E7;
     S_inf_reference = mean(S_xx(idx));
     
    
     
  %% DSH measurement raw data
     S_zz   = (2*pi*f_meas).^2 .* PSD(phase_signal, dt_meas, par.hann_window_switch); % phase_signal =  dphi(t) - dphi(t-tau)   + noise
          
     
  %% maximum likelihood estimation of hidden parameters
   % use hidden time series to find optimal reference values
     
     gaussian = @(x,mu,sigma_sq) 1/sqrt(2*pi*sigma_sq) * exp( -0.5*(x-mu).*(x-mu)/sigma_sq);

   %%%%%%%%%%%%%
   % (1) estimate sigma from hidden S_xixi
         idx = 1E3 <= f_meas & f_meas <= 5E7;

         % data
           f_data = f_meas(idx);
           S_data = S_xixi(idx);
         
         % variable name
           varname = '\sigma';
         
         % proposal function for sigma  
           proposal_sigma = 0.5 * 4E-9;
           proposal_func  = @(x,xprime,sigmasq) gaussian(x,xprime,sigmasq);
  
         % simple uniform prior: enforce positivity
           prior = @(x) prod(double(x>0));

         % number of MCMC steps  
           N_mcmc = 1E5;

         % allocate  
           x    = zeros(N_mcmc, 1);
           logL = zeros(N_mcmc, 1);

         % initialize  
           x(1)    = 1.7E-6;
           logL(1) = loglikelihood_Sxixi(x(1), f_data, S_data);

         % counter  
           accept = 1;
           reject = 0;
  
         % run MCMC  
           for j = 2 : N_mcmc
      
           % update: new proposal      
             xprime = x(j-1) + proposal_sigma * randn(1);

           % evaluate likelihood for new point  
             logL_prime = loglikelihood_Sxixi(xprime, f_data, S_data);

           % accept / reject  
             A = logL_prime + log(prior(xprime)) - logL(j-1) - log(prior(x(j-1)));
             A = min([1,exp(A)]);

             u = rand;

             if u<A % accept
               x(j)    = xprime;
               logL(j) = logL_prime;
               accept  = accept+1;
             else % reject
               x(j)    = x(j-1);
               logL(j) = logL(j-1);
               reject  = reject + 1;
             end

           % plot  
             if mod(j,1E4) == 0 || j == N_mcmc

               figure(601);clf; hold all;
               subplot(1,2,1); hold all;
               h = histogram(x(1:j),'Normalization','pdf','DisplayName','sampled distribution');
                
               x_range = linspace(h.BinEdges(1) - (h.BinEdges(end)-h.BinEdges(1))/2,h.BinEdges(end) + (h.BinEdges(end)-h.BinEdges(1))/2,1001);
               plot(x_range, proposal_func(x_range, xprime, proposal_sigma^2), 'r-','DisplayName','proposal')
               %plot(x_range, prior(x_range), 'g-','DisplayName','prior')

               gauss_x    = 0.5*(h.BinEdges(1:end-1) + h.BinEdges(2:end));
               gauss_mean = trapz(gauss_x, gauss_x.*h.Values);
               gauss_var  = trapz(gauss_x, gauss_x.^2 .* h.Values) - gauss_mean.^2;

               plot(gauss_x, gaussian(gauss_x, gauss_mean, gauss_var), 'b-','DisplayName','dist','LineWidth',2)

               xlabel(varname)
               ylabel('pdf')
                
               sgtitle(['steps: ',num2str(j),...
                        ', acceptance: ',num2str(accept/j),...
                        ', rejection = ',num2str(reject/j)]);
                      
               title(['mean = ',num2str(gauss_mean,'%g'),...
                      ', std = ',num2str(sqrt(gauss_var),'%g')])
                
               ylim([0 max(h.Values)])
               box on
               legend()

               drawnow

             end
              

           end
            

           figure(601); hold all;
           subplot(1,2,2); hold all;
           plot(x, logL, 'ko')
           xlabel(varname)
           ylabel('log likelihood')
           box on
       
         % max likelihood estimate             
           [mu, Sigma] = compute_moments(logL, x);
         
           sigma_ML     = mu;
           sigma_ML_std = sqrt(Sigma);
           
           fprintf(1,'Max likelihood estimate of detector noise\n')
           fprintf(1,'  sigma (mean) = %.4g\n',sigma_ML)
           fprintf(1,'  sigma (std)  = %.4g\n',sigma_ML_std)
                  
                
   %%%%%%%%%%%%%
   % (2) estimate C, nu, S_inf from hidden S_xx
         idx = 1E3 <= f & f <= 5E7;

         % data
           f_data = f(idx);
           S_data = S_xx(idx);
         
         % variable name
           varname = {'C','\nu','S_{\infty}'};
           nvar    = 3;
           
           
         % proposal function for sigma  
           proposal_sigma = 0.1*[1.4E11, 0.044, 2.4];
           proposal_func  = @(x, xprime, sigmasq) gaussian(x, xprime, sigmasq);
  
         % simple uniform prior: enforce positivity
           prior = @(x) prod(double(x>0));

         % number of MCMC steps  
           N_mcmc = 1E5;

         % allocate  
           x    = zeros(N_mcmc, nvar);
           logL = zeros(N_mcmc, 1);

         % initialize  
           x(1,:)  = [1.6e+11 1.51, 468];
           logL(1) = loglikelihood_Sxx(x(1,:), f_data, S_data);

         % counter  
           accept = 1;
           reject = 0;
  
         % run MCMC  
           for j = 2 : N_mcmc
           % update: new proposal      
             xprime = x(j-1,:) + proposal_sigma .* randn(1, nvar);

           % evaluate likelihood for new point  
             logL_prime = loglikelihood_Sxx(xprime, f_data, S_data);

           % accept / reject  
             A = logL_prime + log(prior(xprime)) - logL(j-1) - log(prior(x(j-1,:)));
             A = min([1,exp(A)]);

             u = rand;

             if u<A % accept
               x(j,:)  = xprime;
               logL(j) = logL_prime;
               accept  = accept+1;
             else % reject
               x(j,:)  = x(j-1,:);
               logL(j) = logL(j-1);
               reject  = reject + 1;
             end

           % plot  
             if mod(j,1E4) == 0 || j == N_mcmc

               figure(602);clf; hold all;
               
               for ii = 1 : nvar
               subplot(2,nvar,ii); hold all;
               h = histogram(x(1:j, ii),'Normalization','pdf','DisplayName','sampled distribution');
                
               x_range = linspace(h.BinEdges(1) - (h.BinEdges(end)-h.BinEdges(1))/2,h.BinEdges(end) + (h.BinEdges(end)-h.BinEdges(1))/2,1001);
               plot(x_range, proposal_func(x_range, xprime(ii), proposal_sigma(ii)^2), 'r-','DisplayName','proposal')
               %plot(x_range, prior(x_range), 'g-','DisplayName','prior')

               gauss_x    = 0.5*(h.BinEdges(1:end-1) + h.BinEdges(2:end));
               gauss_mean = trapz(gauss_x, gauss_x.*h.Values);
               gauss_var  = trapz(gauss_x, gauss_x.^2 .* h.Values) - gauss_mean.^2;

               plot(gauss_x, gaussian(gauss_x, gauss_mean, gauss_var), 'b-','DisplayName','dist','LineWidth',2)

               xlabel(varname{ii})
               ylabel('pdf')
                
               sgtitle(['MCMC steps: ',num2str(j),...
                        ', acceptance: ',num2str(accept/j),...
                        ', rejection = ',num2str(reject/j)]);
                      
               title(['mean = ',num2str(gauss_mean,'%g'),...
                      ', std = ',num2str(sqrt(gauss_var),'%g')])
                
               ylim([0 max(h.Values)])
               box on
               legend()

               
               end
               drawnow

               
             end
             
              

           end
             
           
           figure(602); hold all;
           for ii = 1 : nvar            
           subplot(2,nvar,nvar + ii); hold all;
           plot(x(:,ii), logL, 'ko')
           xlabel(varname{ii})
           ylabel('log likelihood')
           box on
           end
       
         % max likelihood estimate      
           [mu, Sigma] = compute_moments(logL, x);
         
           C_ML     = mu(1);
           nu_ML    = mu(2);
           S_inf_ML = mu(3);
           
           C_ML_std     = sqrt(Sigma(1,1));
           nu_ML_std    = sqrt(Sigma(2,2));
           S_inf_ML_std = sqrt(Sigma(3,3));
           
           fprintf(1,'Max likelihood estimate of S_xx signal parameters\n')
           fprintf(1,'  C (mean)     = %.4g\n',C_ML)
           fprintf(1,'  C (std)      = %.4g\n',C_ML_std)
           fprintf(1,'  nu (mean)    = %.4g\n',nu_ML)
           fprintf(1,'  nu (std)     = %.4g\n',nu_ML_std)
           fprintf(1,'  S_inf (mean) = %.4g\n',S_inf_ML)
           fprintf(1,'  S_inf (std)  = %.4g\n',S_inf_ML_std)
       
         
           
%% apply to S_zz
 % estimate C, sigma, nu, S_inf on the fly
 
 
 
 
   %%%%%%%%%%%%%
   
         idx = 1E3 < f_meas & f_meas <= 5E7;
         %idx = 1E3 < f_meas & f_meas <= 1E9;

         % data
           f_data = f_meas(idx);
           S_data = S_zz(idx);
           H_data = Hsq(idx);
           
         % variable name
           varname = {'C','\nu','S_{\infty}','\sigma'};
           nvar    = 4;
                      
         % proposal function for sigma  
           proposal_sigma = 0.1*[7E10, 0.043, 5.8, 7.5E-9];
           proposal_func  = @(x, xprime, sigmasq) gaussian(x, xprime, sigmasq);
  
         % simple uniform prior: enforce positivity
           prior = @(x) prod(double(x>0));

         % number of MCMC steps  
           N_mcmc = 1E5;

         % allocate  
           x    = zeros(N_mcmc, nvar);
           logL = zeros(N_mcmc, 1);

           
           plot_switch= 0;
           
         % initialize  
           x(1,:)  = [1.2E11, 1.45, 450,  1.72E-6];
           logL(1) = loglikelihood_Szz(x(1,:), f_data, S_data, H_data, plot_switch);

         % counter  
           accept = 1;
           reject = 0;
  
         % run MCMC  
           for j = 2 : N_mcmc
               
           % update: new proposal      
             xprime = x(j-1,:) + proposal_sigma .* randn(1, nvar);

           % evaluate likelihood for new point  
             if mod(j,1E4) == 0
                 plot_switch = 1;
             else
                 plot_switch = 0;
             end
             logL_prime = loglikelihood_Szz(xprime, f_data, S_data, H_data, plot_switch);

           % accept / reject  
             A = logL_prime + log(prior(xprime)) - logL(j-1) - log(prior(x(j-1,:)));
             A = min([1,exp(A)]);

             u = rand;

             if u<A % accept
               x(j,:)  = xprime;
               logL(j) = logL_prime;
               accept  = accept+1;
             else % reject
               x(j,:)  = x(j-1,:);
               logL(j) = logL(j-1);
               reject  = reject + 1;
             end

           % plot  
             if mod(j,1E4) == 0 || j == N_mcmc

               figure(702);clf; hold all;
               
               for ii = 1 : nvar
               subplot(2,nvar,ii); hold all;
               h = histogram(x(1:j, ii),'Normalization','pdf','DisplayName','sampled distribution');
                
               x_range = linspace(h.BinEdges(1) - (h.BinEdges(end)-h.BinEdges(1))/2,h.BinEdges(end) + (h.BinEdges(end)-h.BinEdges(1))/2,1001);
               plot(x_range, proposal_func(x_range, xprime(ii), proposal_sigma(ii)^2), 'r-','DisplayName','proposal')
               %plot(x_range, prior(x_range), 'g-','DisplayName','prior')

               gauss_x    = 0.5*(h.BinEdges(1:end-1) + h.BinEdges(2:end));
               gauss_mean = trapz(gauss_x, gauss_x.*h.Values);
               gauss_var  = trapz(gauss_x, gauss_x.^2 .* h.Values) - gauss_mean.^2;

               plot(gauss_x, gaussian(gauss_x, gauss_mean, gauss_var), 'b-','DisplayName','dist','LineWidth',2)

               xlabel(varname{ii})
               ylabel('pdf')
                
               sgtitle(['MCMC steps: ',num2str(j),...
                        ', acceptance: ',num2str(accept/j),...
                        ', rejection = ',num2str(reject/j)]);
                      
               title(['mean = ',num2str(gauss_mean,'%g'),...
                      ', std = ',num2str(sqrt(gauss_var),'%g')])
                
               ylim([0 max(h.Values)])
               box on
               legend()

               
               end
               drawnow

               
             end
             
              

           end
             
           
           figure(702); hold all;
           for ii = 1 : nvar            
           subplot(2,nvar,nvar + ii); hold all;
           plot(x(:,ii), logL, 'ko')
           xlabel(varname{ii})
           ylabel('log likelihood')
           box on
           end
       %%
         % max likelihood estimate      
           [mu, Sigma] = compute_moments(logL, x);
         
           
         % smoothing
           N_smooth = 1000;
           S_xx_smooth = 0;
           f_smooth = 0;
           for i = 1 : N_smooth
           S_xx_smooth = S_xx_smooth + S_xx(i:N_smooth:end);    
           f_smooth    = f_smooth + f(i:N_smooth:end);    
           end
           S_xx_smooth = S_xx_smooth/N_smooth;
           f_smooth    = f_smooth/N_smooth;
           
         % ML results
           S_xx_ML   = mu(1)./(f_data.^mu(2)) + mu(3);
           S_xixi_ML = mu(4)^2 .* f_data.^2;
           S_zz_ML   = H_data.*S_xx_ML + S_xixi_ML;
           
         % PSE reco  
           SNR_ML    = S_xx_ML./S_xixi_ML;
           G_PSE_sq  = SNR_ML./(H_data.*SNR_ML + 1);           
           S_xx_PSE  = G_PSE_sq .* S_data;
           S_xx_inv  = S_data./H_data;
           
           
         % plot  
           figure(703); clf; hold all;
           %
           plot(f(f>0), S_xx(f>0), 'r-','DisplayName','hidden S_{x,x}')
           plot(f_data, S_xx_PSE, 'g-','LineWidth',1,'DisplayName','reco S_{x,x} (PSE)')
           %plot(f_data, S_xx_inv, 'c-','LineWidth',1,'DisplayName','reco S_{x,x} (inv)')

           plot(f_smooth, S_xx_smooth, 'c-','LineWidth',2,'DisplayName','hidden S_{x,x} (smoothed)')
           
           
           
           plot(f_data, S_xx_ML, 'k-','LineWidth',3,'DisplayName','ML est S_{x,x}')
           %}
           %plot(f, movmean(S_xx,5000), 'g-')
           
           %{
           plot(f_data, S_data, 'b-','DisplayName','observed S_{z,z}')
           plot(f_data, S_zz_ML, 'c-','LineWidth',2,'DisplayName','ML est S_{z,z}')
           %}
           
           %{
           plot(f_meas, S_xixi, 'm-','DisplayName','hidden S_{\xi,\xi}')
           plot(f_data, S_xixi_ML, 'k-','LineWidth',3,'DisplayName','ML est S_{\xi,\xi}')
           %}
           set(gca,'XScale','log')
           set(gca,'YScale','log')
           legend()
           
           

           
           
           
  
    
     
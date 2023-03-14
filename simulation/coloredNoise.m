function [F] = coloredNoise(nu, N_OU, N_steps, dt, plot_mode)
   
   % generate a stochastic time series with power-law spectral density
   % S(\omega) \approx \omega^{-\nu} from super position of N_OU 
   % independent Ornstein-Uhlenbeck fluctuators
   %
   % INPUT
   %   nu .......... power law coefficient
   %   N_OU ........ number of fluctuators
   %   N_steps ..... length of time series
   %   dt .......... time step
   %   plot_mode ... switch plotting on (1) or off (0)

%% set regularization
   gamma_0   = 1/(N_steps*dt);
   gamma_inf = 0.1 * 1/dt; % I need to be able to simulate the fastest processes properly, so the rate must satisfy gamma_inf*dt < 1
   
   
   
%% generate rates  gamma_j
   % check
     if nu >= 2
       error('nu too large')
     elseif nu <= 0
       error('nu too small')
     end

     
   % inverse transform sampling  
     U = linspace(0, 1, N_OU);
     
     if nu == 1
       % special case nu = 1
         C_nu  = 1/log(gamma_inf/gamma_0);         
         gamma = gamma_0 * exp(U/C_nu);
     else
       % for 0 < nu < 2 (except for nu=1)    
         C_nu  = (1-nu)/(gamma_inf^(1-nu) - gamma_0^(1-nu));
         gamma = ((1-nu)/C_nu * U + gamma_0^(1-nu)).^(1/(1-nu));   
     end
     
     
     %figure(2413);clf;hold all;
     %plot(1:N_OU,gamma,'ko')
     %set(gca,'YScale','log')


     
%% set scaling parameter to normalize PSD: S_FF = 1/omega^nu
   scaling_A = sin(pi*nu/2)/(C_nu * pi);
     
   
%% simulate Ornstein-Uhlenbeck flucturators  
   X = zeros(N_OU,N_steps);
   for i = 1 : N_steps-1
     X(:,i+1) = (1 - dt*gamma').*X(:,i) + sqrt(2*dt*gamma').*randn(N_OU,1);
   end
     
 % generate colored noise  
   F = sqrt(scaling_A) * sum(X,1)/sqrt(N_OU);
     
   
%% plots    
   if plot_mode == 1     
     
     figure(101);clf;hold all;
   
       % plot histogram of rates
       subplot(2,2,1); hold all;
          histogram(gamma,'Normalization','pdf')
          range = linspace(gamma_0,gamma_inf,10001);
          dist  = @(gamma) C_nu * gamma.^(-nu) .* heaviside(gamma-gamma_0).*heaviside(gamma_inf-gamma);
          plot(range, dist(range), 'r-','LineWidth',2)
          title('distribution of rates')
          xlabel('\gamma')
          ylabel('\rho(\gamma)')
          %set(gca,'XScale','log')          
          set(gca,'YScale','log')
          box on;
          
          trapz(range,dist(range))
     
       % OU fluctuators
       subplot(2,2,2); hold all;
         for i = 1 : N_OU
           plot([0:N_steps-1]*dt, 2*i + X(i,:),'-','LineWidth',0.5)
         end     
         plot([0:N_steps-1]*dt, F,'k-','LineWidth',2)
         box on;
         xlabel('time t [s]')
         ylabel('noise F [1/s]')
         title('time series')
         xlim([0, (N_steps-1)*dt])     
         
       % ACF
       subplot(2,2,3);hold all;
         maxlag = round(N_steps/32);
         [C_FF, lags] = xcorr(F,maxlag,'normalized');
         tau  = lags * dt;
         plot(tau, C_FF,'k-','LineWidth',2,'DisplayName','numeric')
         box on
         xlabel('\tau')
         ylabel('C_{F,F}(|\tau|)')
         title('autocorrelation function (normalized)')
         xlim([min(tau),max(tau)])        
         
         % compute and plot analytic autocorrelation function
           C_FF_analytic = zeros(1,length(tau));
           for i = 1 : N_OU
           C_FF_analytic = C_FF_analytic + exp(-gamma(i)*abs(tau));
           end
           C_FF_analytic = scaling_A * C_FF_analytic/N_OU;             % normalization           
           C_FF_analytic = C_FF_analytic/C_FF_analytic(find(lags==0)); % normalize to 1 at tau=0
           plot(tau,C_FF_analytic,'r--','LineWidth',2,'DisplayName','analytic')
         
         legend()
         
       % PSD
       subplot(2,2,4);hold all;
         f             = [-N_steps/2 : N_steps/2-1]*1/(N_steps*dt);
         S_FF          = dt/N_steps * abs(fftshift(fft(F))).^2;
         %S_FF_analytic = C_nu*pi/sin(pi*nu/2) *  (2*pi*f).^(-nu);
         S_FF_analytic = (2*pi*f).^(-nu);
         
         plot(f, S_FF,'k-','LineWidth',1,'DisplayName','numeric')
         plot(f, S_FF_analytic,'r-','LineWidth',2,'DisplayName','analytic')        

         plot([1 1]*gamma_inf, [min(S_FF) max(S_FF)],'k--','LineWidth',1,'DisplayName','max rate \gamma_\infty')
         
         % linear least squares fit
         %{
           range = (N_steps/2+2) : N_steps;
           [p,S] = polyfit(log(2*pi*f(range)),log(S_FF(range)),1);
           C_fit  = exp(p(2));
           nu_fit = -p(1);
           plot(f(range), C_fit./(2*pi*f(range).^nu_fit),'g--','DisplayName','linear fit')
         %}
         
         % nonlinear least squares fit  
         %{
           p   = optimvar('p',2)
           fun = p(1)/(2*pi*f(range)^p(2))
           obj = sum((fun - S_FF(range)).^2);
           
           lsqproblem = optimproblem("Objective",obj);

           x0.p = [1, nu]

           show(lsqproblem)

           [sol,fval] = solve(lsqproblem,x0)
         %}
         set(gca,'XScale','log')
         set(gca,'YScale','log')
         xlabel('frequency \it{f}')
         ylabel('\it{S_{F,F}(f)}')
         box on
         axis tight    
         legend()
         
         title('power spectral density')
         
         
         sgtitle(['colored noise: mean(F) = ',num2str(mean(F)),', var(F) = ',num2str(var(F))])
       


       
       
   end
   
   
   
   
   
  
   
   return;
   %% dFdF
   
   
     dF = [diff(F)];
     
     mean_dFdF = mean(dF.*dF)
   
     
     %{
     mean_gamma1 = mean(gamma)
     
     rho = @(gamma) gamma.^(-nu);
     gamma_range = logspace(log10(gamma_0), log10(gamma_inf), 10001);
     mean_gamma2 = trapz(gamma_range,gamma_range.*rho(gamma_range))     
     mean_gamma3 = C_nu/(2-nu) * (gamma_inf^(2-nu) - gamma_0^(2-nu))
     
     
     mean_dFdF_trial1 = 2*scaling_A*mean_gamma1 * dt
     mean_dFdF_trial2 = 2*scaling_A*mean_gamma2 * dt
     mean_dFdF_trial3 = 2*scaling_A*mean_gamma3 * dt
     
     %}
     
     sum(gamma)/N_OU
     
     mean(gamma)
     
     scaling_A
     
     
     mean_dFdF_trial = scaling_A * 2 * sum(gamma)/N_OU * dt
     
     
  error_quotient =mean_dFdF/mean_dFdF_trial
  
  
  
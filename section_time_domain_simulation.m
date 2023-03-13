   simulate_nonlinear_system     = 1;
   simulate_linear_system        = 0;

   plot_time_series              = 1;
   plot_power_spectral_densities = 1;
  
   
   dt  = 50E-12;         
   
   % discrete delay time
   l   = round(par.tau_d/dt);
   
   N_steps = 100*l + l; % add one interferometer delay period (will be deleted later on when simulating the experiment)
   
   
   
         
   t = [0:N_steps-1]*dt;
   f = [-N_steps/2 : N_steps/2-1]/(N_steps*dt);
  
   
   % print information
     fprintf(1,'\n')
     fprintf(1,'=== time domain simulation ===\n')
     fprintf(1,'    num steps = %d (%.4e)\n',N_steps,N_steps)
     fprintf(1,'           dt = %.4e \n',dt)
     fprintf(1,'            T = %.4e \n',N_steps*dt)
     fprintf(1,'        f_min = %.4e \n',1/(N_steps*dt))
     fprintf(1,'        f_max = %.4e \n',1/(2*dt))
     fprintf(1,'\n')
   
   % select pump current
     par.I        = I_range(end);            % pick the last one from the PSD section to reuse the exact PSDs in the plot below
     [P_ss, N_ss] = steady_state(par);
     
   % CW frequency  
     Omega_CW     = 0.5*par.alpha_H* par.Gamma*par.vg*func_g(par,P_ss,N_ss);

  %% generate white noise (twelve processes)
     fprintf(1,'    generate white noise ............... ');
     tic
     dW = sqrt(dt) * randn(12,N_steps);
     toc

  %% generate colored noise
     fprintf(1,'    generate colored noise ............. ');
     tic
     F      = zeros(3, N_steps); %F = [F_P;F_phi;F_N];
     F(1,:) = coloredNoise(par.nu(1), par.N_F(1), N_steps, dt, 0);
     F(2,:) = coloredNoise(par.nu(2), par.N_F(2), N_steps, dt, 0);
     F(3,:) = coloredNoise(par.nu(3), par.N_F(3), N_steps, dt, 0);
     toc
     
     
  %% nonlinear system  
     if simulate_nonlinear_system == 1
     fprintf(1,'    simulate nonlinear system .......... ');
     X0 = [P_ss, 0, N_ss]; % start on steady state
     tic
     sol_nlin = sim_nlin(par, dt, N_steps, X0, dW, F);
     toc
     end
     
  %% linearized system
     if simulate_linear_system == 1
     fprintf(1,'    simulate linearized system ......... ');
     X0 = [0, 0, 0]; % start on steady state
     tic
     sol_lin = sim_lin(par, dt, N_steps, X0, dW, F);
     toc
     end
     
  %% plot time series
     if plot_time_series == 1
     fprintf(1,'    plot time series ................... ');     
     tic
     
     figure(10); clf; hold all;
     
       sgtitle(['nonlinear vs. linearized dynamics (time series) | I = ',num2str(par.I*1E3),' mA'])  
       
       subplot(2,2,1); hold all;
         if simulate_nonlinear_system == 1            
         plot(t, sol_nlin.P - P_ss, 'k-','LineWidth',2,'DisplayName','\deltaP(t) (nonlinear)')
         end
         if simulate_linear_system == 1
         plot(t, sol_lin.P,         'r-.','LineWidth',1,'DisplayName','\deltaP(t) (linearized)')
         end
         box on;
         xlabel('t [s]')
         ylabel('fluctuations \deltaP')
         axis tight
         legend('Location','southwest')
         title('photon number fluctuations')
       
       subplot(2,2,2); hold all;
         if simulate_nonlinear_system == 1
         plot(t, sol_nlin.phi - Omega_CW*t, 'k-','LineWidth',2,'DisplayName','\delta\phi(t) (nonlinear)')
         end
         if simulate_linear_system == 1
         plot(t, sol_lin.phi,               'g-.','LineWidth',1,'DisplayName','\delta\phi(t) (linearized)')
         end
         xlabel('t [s]')
         ylabel('fluctuations \delta\phi')
         box on;
         axis tight
         legend('Location','southwest')
         title('phase fluctuations')
                            
       subplot(2,2,3); hold all;
         if simulate_nonlinear_system == 1
         plot(t, sol_nlin.N - N_ss, 'k-','LineWidth',2,'DisplayName','\deltaN(t) (nonlinear)')
         end
         if simulate_linear_system == 1
         plot(t, sol_lin.N,         'b-.','LineWidth',2,'DisplayName','\deltaN(t) (linearized)')
         end
         box on;
         xlabel('t [s]')
         ylabel('fluctuations \deltaN')
         axis tight
         legend('Location','southwest')  
         title('carrier number fluctuations')
         
       subplot(2,2,4); hold all;
         if simulate_nonlinear_system == 1
         plot(t, [0,diff(sol_nlin.phi - Omega_CW*t)/dt], 'k-','LineWidth',2,'DisplayName','\delta\omega(t) (nonlinear)')         
         end
         if simulate_linear_system == 1
         plot(t, [0,diff(sol_lin.phi)/dt],                'c-.','LineWidth',1,'DisplayName','\delta\omega(t) (linearized)')
         plot(t, movmean([0,diff(sol_lin.phi)/dt], 1000), 'k-','LineWidth',2,'DisplayName','\delta\omega(t) (moving average)')
         end
         xlabel('t [s]')
         ylabel('frequency fluctuations')
         box on;
         axis tight
         legend('Location','southwest')    
         title('frequency fluctuations')
         xlim([0,1E-8])
         
         
         drawnow       
       
     toc  
     end
       

     
  %% plot power spectral densities
  
     if plot_power_spectral_densities == 1
     fprintf(1,'    plot power spectral densities ...... ');     
     
     tic
     
     figure(11); clf; hold all;
     
       sgtitle(['nonlinear vs. linearized dynamics (power spectra) | I = ',num2str(par.I*1E3),' mA'])  
       
       subplot(3,1,1); hold all;
         if simulate_nonlinear_system == 1
         plot(f, PSD(sol_nlin.P - P_ss, dt, par.hann_window_switch)/P_ss^2, 'k-','LineWidth',1,'DisplayName','RIN(f) (nonlinear)')
         end
         if simulate_linear_system == 1
         plot(f, PSD(sol_lin.P, dt, par.hann_window_switch)/P_ss^2,         'r-.','LineWidth',1,'DisplayName','RIN(f) (linearized)')
         end         
         plot(f_analytic,S_dPdP/P_ss^2,                                     'c-','LineWidth',2,'DisplayName','RIN(f) (exact)')
         box on;
         xlabel('f [Hz]')
         ylabel('S_{\deltaP,\deltaP}/P^2')
         axis tight
         set(gca,'XScale','log')
         set(gca,'YScale','log')
         legend('Location','southwest')
         
       
       subplot(3,1,2); hold all;
         if simulate_nonlinear_system == 1
         plot(f, (2*pi*f).^2 .* PSD(sol_nlin.phi - Omega_CW*t, dt, par.hann_window_switch), 'k-','LineWidth',1,'DisplayName','FN PSD (nonlinear)')
         end
         if simulate_linear_system == 1
         plot(f, (2*pi*f).^2 .* PSD(sol_lin.phi, dt, par.hann_window_switch),               'g-','LineWidth',1,'DisplayName','FN PSD (linearized)')
         end
         plot(f_analytic,S_domegadomega,                                                    'c-','LineWidth',2,'DisplayName','FN PSD(f) (exact)')
         xlabel('f [Hz]')
         ylabel('S_{\deltaf,\deltaf}')
         box on;
         axis tight
         set(gca,'XScale','log')
         set(gca,'YScale','log')
         legend('Location','southwest')
         
              
       subplot(3,1,3); hold all;
         if simulate_nonlinear_system == 1
         plot(f, PSD(sol_nlin.N - N_ss, dt, par.hann_window_switch)/N_ss^2, 'k-','LineWidth',1,'DisplayName','rel. carrier noise PSD (nonlinear)')
         end
         if simulate_linear_system == 1
         plot(f, PSD(sol_lin.N, dt, par.hann_window_switch)/N_ss^2,         'b-','LineWidth',1,'DisplayName','rel. carrier noise PSD (linearized)')
         end
         plot(f_analytic,S_dNdN/N_ss^2,                                     'c-','LineWidth',2,'DisplayName','rel. carrier noise PSD (exact)')
         box on;
         xlabel('f [Hz]')
         ylabel('S_{\deltaN,\deltaN}/N^2')
         axis tight
         set(gca,'XScale','log')
         set(gca,'YScale','log')         
         legend('Location','southwest') 
         
      
         
     toc
     end
     
     
     
     
     
         
       

     

    %% plot for Fig. 4 in Opt. Express paper
         
         
       figure(100); clf; hold all;

         subplot(1,2,1); hold all;
       
         if simulate_nonlinear_system == 1         
         aux = plot(f, (2*pi*f).^2 .* PSD(sol_nlin.phi - Omega_CW*t, dt, par.hann_window_switch), '-','Color',[1 1 1]*0.6,'LineWidth',1,'DisplayName','FN PSD (nonlinear)');
         uistack(aux,'bottom')
         end
         if simulate_linear_system == 1
         aux = plot(f, (2*pi*f).^2 .* PSD(sol_lin.phi, dt, par.hann_window_switch),               '-','Color',[1 1 1]*0.6,'LineWidth',1,'DisplayName','FN PSD (linearized)');
         uistack(aux,'bottom')
         end
         
         plot(f_analytic, S_domegadomega,        'k-','LineWidth',3,'DisplayName','white + colored noise')
         plot(f_analytic, S_domegadomega_WN,     'c--','LineWidth',2,'DisplayName','white noise')
         plot(f_analytic, S_domegadomega_CN_P,   'r:','LineWidth',2,'DisplayName','colored noise P')
         plot(f_analytic, S_domegadomega_CN_phi, 'g:','LineWidth',2,'DisplayName','colored noise \phi')
         plot(f_analytic, S_domegadomega_CN_N,   'b:','LineWidth',2,'DisplayName','colored noise N')
         legend('Location','southwest','AutoUpdate','off')
                  
         box on
         axis tight
         set(gca,'XScale','log')
         set(gca,'YScale','log')
         xlabel('f [Hz]')
         ylabel('S_{\delta\omega, \delta\omega} [Hz^2/Hz]')
         title('frequency noise')
         
         xlim([1E3, 5E9])
         ylim([1E-6, 1E8])
         
         xzoom_min = 1E7;
         xzoom_max = 5E9;
         
         yzoom_min = 800;
         yzoom_max =  10000;
         
         plot([xzoom_min xzoom_max xzoom_max xzoom_min xzoom_min], [yzoom_min yzoom_min yzoom_max yzoom_max yzoom_min],'m-')
         
       % zoom         
         subplot(1,2,2); hold all;
         
         skip = 25;
       
         idx_range = find(f==0) : skip : length(f);
         
         if simulate_nonlinear_system == 1    
         S = (2*pi*f).^2 .* PSD(sol_nlin.phi - Omega_CW*t, dt, par.hann_window_switch);
         aux = plot(f(idx_range), S(idx_range), '-','Color',[1 1 1]*0.6,'LineWidth',1,'DisplayName','FN PSD (nonlinear)');
         uistack(aux,'bottom')
         end
         if simulate_linear_system == 1
         S = (2*pi*f).^2 .* PSD(sol_lin.phi, dt, par.hann_window_switch);
         aux = plot(f(idx_range), S(idx_range),               '-','Color',[1 1 1]*0.6,'LineWidth',1,'DisplayName','FN PSD (linearized)');
         uistack(aux,'bottom')
         end
         
         plot(f_analytic, S_domegadomega,        'k-','LineWidth',3,'DisplayName','white + colored noise')
         plot(f_analytic, S_domegadomega_WN,     'c--','LineWidth',2,'DisplayName','white noise')
         plot(f_analytic, S_domegadomega_CN_P,   'r:','LineWidth',2,'DisplayName','colored noise P')
         plot(f_analytic, S_domegadomega_CN_phi, 'g:','LineWidth',2,'DisplayName','colored noise \phi')
         plot(f_analytic, S_domegadomega_CN_N,   'b:','LineWidth',2,'DisplayName','colored noise N')
                  
         box on
         axis tight
         set(gca,'XScale','log')
         set(gca,'YScale','log')
         xlabel('f [Hz]')
         ylabel('S_{\delta\omega, \delta\omega} [Hz^2/Hz]')
         title('frequency noise')
         
         xlim([xzoom_min, xzoom_max])
         ylim([yzoom_min, yzoom_max])
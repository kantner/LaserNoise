fprintf(1,'\n')
fprintf(1,'=== simulate SDH measurement ===\n')

   % shift tau_d on time grid  
     l            = round(par.tau_d/dt);
     par.tau_d    = l*dt;                 
     
     
  %% generate I-Q data     
     % choose which simulation data shall be used to simulate the measurement
  
     time_series_model = 0;
     % 0 = nonlinear model
     % 1 = linearized model
     
     fprintf(1,'    generate I-Q data .................. ')
     tic
     switch(time_series_model)
       case 0  % generate I-Q data from sol_nlin
         [I, Q, I_no_noise, Q_no_noise] = generate_IQdata(sol_nlin.P      , sol_nlin.phi            , l, dt, par.deltaOmega, par.eta_det, par.eta_out, par.sigma_IQ);
       case 1  % generate I-Q data from sol_lin
         [I, Q, I_no_noise, Q_no_noise] = generate_IQdata(P_ss + sol_lin.P, Omega_CW*t + sol_lin.phi, l, dt, par.deltaOmega, par.eta_det, par.eta_out, par.sigma_IQ);
     end
     toc
          
     % time and frequency for measurement data (is shorter since leading delay period is erased): N_meas = N_steps - l
     N_meas = length(I);    
     t_meas = [0:N_meas-1]*dt;
     f_meas = [-N_meas/2 : N_meas/2-1]/(dt*N_meas);

     dt_meas = dt;
     l_meas  = l;              
     
     
  %% filter out high-frequency components
  %{
   % coarse graining of time series
     f_max_target = 5E8;
     n_coarse = round(1/(2*dt*f_max_target));
     n_coarse = 2^round(log2(n_coarse));
     
     coarse_graining_switch = 0; % 0 = off (no coarse graining) , 1 = on (reduce sampling rate) 
     
     switch(coarse_graining_switch)
       case 0
         l_meas  = l;         
         dt_meas = t_meas(2) - t_meas(1);
         
       case 1
         mode_coarse = 2;
         % 1: fft
         % 2: averaging
         % 3: stroboscopic
       
       % overwrite 
         [I,          t_meas_coarse] = reduce_sample_rate(I,          t_meas, n_coarse, mode_coarse);
         [Q,          t_meas_coarse] = reduce_sample_rate(Q,          t_meas, n_coarse, mode_coarse);     
         [I_no_noise, t_meas_coarse] = reduce_sample_rate(I_no_noise, t_meas, n_coarse, mode_coarse);
         [Q_no_noise, t_meas_coarse] = reduce_sample_rate(Q_no_noise, t_meas, n_coarse, mode_coarse);
                 
         N_meas     = length(I);
         t_meas     = t_meas_coarse;
         dt_meas    = mean(diff(t_meas_coarse));
         f_meas     = [-N_meas/2 : N_meas/2-1]/(dt_meas*N_meas);
     
         l_meas     = round(par.tau_d/dt_meas);
       otherwise  
         error('not implemented');  
     end
  %}
     
  %%   
     
   figure(20); clf; hold all;
       
   sgtitle('Measurement: Time Series')

   subplot(3,1,1); hold all;
     plot(t_meas, I,          'r-', 'LineWidth',2,'DisplayName','I (with meas. noise)')
     plot(t_meas, Q,          'b-', 'LineWidth',2,'DisplayName','Q (with meas. noise)')
     plot(t_meas, I_no_noise, 'r-.','LineWidth',1,'DisplayName','I (without meas. noise)')
     plot(t_meas, Q_no_noise, 'b-.','LineWidth',1,'DisplayName','I (without meas. noise)')
     title('I and Q data')
     legend()
     xlabel('time [s]')
     ylabel('I(t), Q(t) [W]')
     box on;
     axis tight
   
     
     
     
   % extract RIN-type signal
     power_norm            = mean(I.^2 + Q.^2);             % approximate (P_ss * par.eta_det * par.eta_out)^2
     power_signal          = (I.^2 + Q.^2)/power_norm - 1;  % approximate [deltaP(t) + deltaP(t-tau)]/P0
     
     power_norm_no_noise   = mean(I_no_noise.^2 + Q_no_noise.^2);
     power_signal_no_noise = (I_no_noise.^2 + Q_no_noise.^2)/power_norm_no_noise - 1;
     
   % extract exact power noise
     power_noise           = power_signal - power_signal_no_noise;
     

   subplot(3,1,2); hold all;  
     plot(t_meas, power_signal         ,'r-' ,'LineWidth',2,'DisplayName','power signal (with meas. noise)')
     plot(t_meas, power_signal_no_noise,'r-.','LineWidth',1,'DisplayName','power signal (without meas. noise)')
     xlabel('time [s]')
     ylabel('I^2 + Q^2 [(W)^2]')
     box on;
     title('sum of power fluctuations')
     axis tight
     legend();
   
   

   % extract phase signal (for FN reconstruction)
     phase_signal          = atan(Q./I)                   + par.deltaOmega * t_meas - Omega_CW * par.tau_d;   % approximate \delta\phi(t) - \delta\phi(t-tau)
     phase_signal_no_noise = atan(Q_no_noise./I_no_noise) + par.deltaOmega * t_meas - Omega_CW * par.tau_d; 
     
          
   % correct for phase jumps modulo [-pi,pi]  
     fprintf(1,'    correct phase jumps (with noise) ... ')
     tic
     phase_signal          = correct_phase_jumps(phase_signal);
     toc
     
     fprintf(1,'    correct phase jumps (no noise) ..... ')
     tic
     phase_signal_no_noise = correct_phase_jumps(phase_signal_no_noise);
     toc
     
     
   % extract exact phase noise
     phase_noise = phase_signal - phase_signal_no_noise;

   subplot(3,1,3); hold all;      
     plot(t_meas, phase_signal          ,'g-' ,'LineWidth',2,'DisplayName','phase signal (with meas. noise)')
     plot(t_meas, phase_signal_no_noise ,'g-.','LineWidth',1,'DisplayName','phase signal (without meas. noise)')
     xlabel('time [s]')
     ylabel('arctan(Q/I) + \Delta\omegat - \Omega \tau_d')
     box on;
     title('difference of phase fluctuations')
     axis tight
     legend();
   
   % plot measurement results   
     fprintf(1,'    plot power spectral densities ...... ')
     tic
     
     
   figure(21); clf; hold all;
   
   sgtitle('Measurement: Power Spectral Densities')
     
     % theory     
       S_powerNoise = (2*par.sigma_IQ*dt)^2/power_norm;
     % empirical
     %  S_powerNoise = mean(PSD(power_noise, dt_meas));
   
     subplot(2,1,1); hold all;
     plot(f_meas, PSD(power_signal, dt_meas)           ,'r-', 'LineWidth',1,'DisplayName','power signal (with meas. noise)')
     plot(f_meas, PSD(power_signal_no_noise, dt_meas)  ,'k-', 'LineWidth',1,'DisplayName','power signal (without meas. noise)')
     plot(f_meas, PSD(power_noise, dt_meas)            ,'m-', 'LineWidth',1,'DisplayName','meas. noise')
     plot(f_meas, S_powerNoise * ones(size(f_meas))    ,'m-.','LineWidth',4,'DisplayName','meas. noise (analytic)')
     set(gca,'XScale','log')
     set(gca,'YScale','log')
     xlabel('f [Hz]')
     ylabel('PSD (power signal)')
     axis tight
     box on
     legend();
     title('PSD power signal')
     
     % theory
       S_phaseNoise = 0.25 * S_powerNoise;
     % empirical
     %  S_phaseNoise = mean(PSD(phase_noise, dt_meas));
     
     subplot(2,1,2); hold all;
     plot(f_meas, PSD(phase_signal, dt_meas)            ,'g-', 'LineWidth',1,'DisplayName','phase signal (with meas. noise)')
     plot(f_meas, PSD(phase_signal_no_noise, dt_meas)   ,'k-', 'LineWidth',1,'DisplayName','phase signal (without meas. noise)')
     plot(f_meas, PSD(phase_noise, dt_meas)             ,'m-', 'LineWidth',1,'DisplayName','meas. noise')
     plot(f_meas, S_phaseNoise * ones(size(f_meas))     ,'m-.','LineWidth',4,'DisplayName','meas. noise (analytic)')
     set(gca,'XScale','log')
     set(gca,'YScale','log')
     xlabel('f [Hz]')
     ylabel('PSD (phase signal)')
     axis tight
     box on
     legend();
     title('PSD phase signal')
     
     toc
     
     
     
     
  fprintf(1,'\n') 
  fprintf(1,'=== filtering of SDH signal ===\n')
  % The reconstruction is for the FN-PSD only.
  
  %% convolution function   
   % convolution function h(t) and Fourier transform H(\omega) for phase signal       
     h                   = zeros(N_meas,1);
     h(1)                = 1;  % instantaneous signal
     h(N_meas+1-l_meas)  = -1; % retarded signal
     
     H                   = transpose(fftshift(fft(h)));
     
     %Hsq = abs(H).^2;
     
     % direct
     Hsq = 2*(1 - cos(2*pi*f_meas.*par.tau_d));
     
     
     
        
  %% hidden signal and noise
     S_xx   = (2*pi*f).^2 .* PSD(sol_nlin.phi - Omega_CW*t, dt, par.hann_window_switch);
     S_xixi = (2*pi*f_meas).^2 .* PSD(phase_noise, dt_meas, par.hann_window_switch);     
  
 
  %% DSH measurement raw data
     S_zz   = (2*pi*f_meas).^2 .* PSD(phase_signal, dt_meas, par.hann_window_switch); % phase_signal =  dphi(t) - dphi(t-tau)   + noise
                    
        
  %% inverse filter
     S_xhatxhat_inverse = S_zz./Hsq;  
     
  %% estimate parameters: simple fitting
  fprintf(1,'    estimate parameters ....')
  tic

  
   % assume functional form
   %         C
   % S_xx = ---- + S_inf      S_xixi = sigma^2 * f^2
   %        f^nu
   
   figure(32); clf; hold all;
     sgtitle('parameter estimation')
   
   %%%%%%%%%%%%%%%%%%%%%%%
   % (1) estimate sigma from values of S_zz at pole frequencies
   
         subplot(1,2,1); hold all;
         box on;
         xlabel('f (Hz)')
         ylabel('S_{z,z} at poles')
      
       % find pole indices        
         idx = find(abs(H).^2 < 1E-12 & f_meas > 0);
        
       % plot values of S_zz at pole indices          
         plot(f_meas(idx), S_zz(idx), 'bo', 'DisplayName', 'S_{z,z} (f_{pole})')
         
       % quadratic least squares fit  
         X = f_meas(idx);
         Y = S_zz(idx);
         sigma_fit = sqrt( sum(X.^2 .* Y)./sum(X.^4) );
         
         plot(X, sigma_fit^2 .* X.^2, 'm-','LineWidth',2, 'DisplayName', 'quadratic fit: \sigma^2 f^2')
         
         legend('Location','northwest')
         title(['fit \sigma = ',num2str(sigma_fit)])
         
         %{
         fprintf(1, 'detector noise level estimation\n')
         fprintf(1, '  sigma = %.4g\n\n',sigma_fit)
         %}
         
   % (2) estimate C and nu from values of S_xhathat_inverse in low frequency regime
         subplot(1,2,2); hold all;
         box on;
         xlabel('f (Hz)')
         ylabel('S_{x,x}^{rec} (inverse)')
         legend();
         
         %idx = find( 5E3 < f_meas & f_meas < 8E4 ), ...
       % select frequency range for fitting
         idx = [find( 1E3 <= f_meas & f_meas <= 9E4 ),...
                find( 1.1E5 <= f_meas & f_meas <= 1.9E5),...
                %find( 2.1E5 <= f_meas & f_meas <= 2.9E5)...
               ];
         
         
         plot(f_meas(f_meas>0), S_xhatxhat_inverse(f_meas>0), 'r-', 'DisplayName','S_{x,x}^{rec} (inv)')
         plot(f_meas(idx), S_xhatxhat_inverse(idx), 'b-', 'DisplayName','S_{x,x}^{rec} (inv): fitting range')
         set(gca,'XScale','log')
         set(gca,'YScale','log')
         
         plot(f_meas(f_meas>0), sigma_fit^2 * f_meas(f_meas>0).^2, 'm-','LineWidth',2, 'DisplayName','noise fit')
         
     
       % linear least squares  
         xdata = log10(f_meas(idx));
         ydata = log10(S_xhatxhat_inverse(idx));
   
         A = [mean(xdata.^2), mean(xdata);
              mean(xdata),    1];
         b = [mean(xdata.*ydata);
              mean(ydata)];
   
         parameters = A\b; % solve least squares problem
   
         C_linLSQ  = 10^parameters(2);
         nu_linLSQ = -parameters(1);
     
         %{
         fprintf(1,'  linear least fit ...\n')
         fprintf(1,'    C  = %.4e\n',C_linLSQ)
         fprintf(1,'    nu = %.4f\n',nu_linLSQ)
         %}
         
         plot(f_meas(f_meas>0), C_linLSQ./( f_meas(f_meas>0).^nu_linLSQ) ,'g-', 'DisplayName','linear LSQ','LineWidth',2)
         
     
       % nonlinear least squares    
         xdata = f_meas(idx);
         ydata = S_xhatxhat_inverse(idx);
      
         func = @(p, xdata) 10^p(2) .* xdata.^p(1); % fit function
   
         p0 = [parameters(1), parameters(2)]; % use linear least squares result as starting guess
   
         opts = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
         lower_bound = [];
         upper_bound = [];
         p = lsqcurvefit(func, p0, xdata, ydata, lower_bound, upper_bound, opts);
   
         C_nlinLSQ   = 10^p(2);
         nu_nlinLSQ = -p(1);
     
         %{
         fprintf(1,'  nonlinear least fit ...\n')
         fprintf(1,'    C  = %.4e\n',C_nlinLSQ)
         fprintf(1,'    nu = %.4f\n',nu_nlinLSQ)   
         %}
         
         plot(f_meas(f_meas>0), C_nlinLSQ ./ ( f_meas(f_meas>0).^nu_nlinLSQ ),'b-', 'DisplayName','nonlinear LSQ','LineWidth',2)
         ylim([1E0 1E8])
     
         
         toc
         
         
  %% PSE filter: Minimization
     fprintf(1,'    minimize objective function ... ')
     tic

  
   % Sweep over S_inf
     S_inf_range       = logspace(1,6,401);
     % this is now set outside the loop
     
     
     consistency_error = zeros(1,length(S_inf_range));
     
   % set range for error evaluation 
     idx = (1E3 < f_meas) & (f_meas <= 5E7); % for larger frequencies, the footprint of the RO peak becomes visible
     
     C_est     = C_nlinLSQ;
     nu_est    = nu_nlinLSQ;
     sigma_est = sigma_fit;
     
     plot_switch = 0;
     
     S_xx_model     = @(S_inf) C_est./f_meas(idx).^nu_est + S_inf;
     SNR            = @(S_inf) S_xx_model(S_inf)./(sigma_est^2 * f_meas(idx).^2);
     G_PSE_sq       = @(S_inf) SNR(S_inf)./( Hsq(idx).*SNR(S_inf) + 1 );
     S_xhatxhat_PSE = @(S_inf) G_PSE_sq(S_inf) .* S_zz(idx);  
     objective_func = @(S_inf) mean( (S_xhatxhat_PSE(S_inf) - S_xx_model(S_inf))./S_xx_model(S_inf) );
     
     for i = 1 : length(S_inf_range)
                
       consistency_error(i) = objective_func(S_inf_range(i));
       
       if plot_switch == 1
       figure(200);clf; hold all;
       plot(f(f>0), S_xx(f>0), 'r-')
       plot(f_meas(idx), S_xhatxhat_PSE(S_inf_range(i)), 'g-')
       plot(f_meas(idx), S_xx_model(S_inf_range(i)), 'k-','LineWidth',2)       
       set(gca,'XScale','log')
       set(gca,'YScale','log')       
       end
       
       

       
     end
    
   % find minimum via bisection
     opts      = optimset('TolX',1E-36,'Display','none');
     S_inf_est = fzero(objective_func, [1E1 1E6], opts);
     
   
     [S_inf_range, order] = sort([S_inf_range, S_inf_est]);
     consistency_error    = [consistency_error, objective_func(S_inf_est)];
     consistency_error    = consistency_error(order);
     
     %
     figure(201);clf;hold all;
       plot(S_inf_range, consistency_error.^2, 'r-','DisplayName','consistency error')       
       set(gca,'XScale','log')
       set(gca,'YScale','log')
       title('Minimization of consistency error')
       ylabel('objective function D(S_{\infty}')
       xlabel('S_{\infty}')
       box on
       legend('Location','southeast')
       ylim([5E-7 5E0])

       drawnow
       
    toc
%% plot comparison of all filters at optimized parameters
   fprintf(1,'    plotting ...................... ')
   tic


   S_xx_est       = C_est./f_meas.^nu_est + S_inf_est;
   S_xixi_est     = (sigma_est^2 * f_meas.^2);
   SNR            = S_xx_est./S_xixi_est;
   G_PSE_sq       = SNR./( Hsq.*SNR + 1);
   G_Wiener_sq    = Hsq.*SNR.^2 ./( Hsq.*SNR + 1).^2;
   
   S_zz_est       = Hsq.*S_xx_est + S_xixi_est;
       
   S_xhatxhat_PSE    = G_PSE_sq    .* S_zz;
   S_xhatxhat_Wiener = G_Wiener_sq .* S_zz;
   
   figure(300); clf; hold all;
   
     xmin = 1E3;
     xmax = 1E8;
     ymin = 5E-3;
     ymax = 5E6;

     
     
     subplot(2,2,1); hold all     
     plot(f(f>0)          , S_xx(f>0), 'r-','LineWidth',1,'DisplayName','S_{x,x} (hidden)')
     plot(f_meas(f_meas>0), S_xixi(f_meas>0), 'm-','LineWidth',1,'DisplayName','S_{\xi,\xi} (hidden)')     
     plot(f_meas(f_meas>0), S_zz(f_meas>0), 'b-','LineWidth',1,'DisplayName','S_{z,z}')
     
     plot(f_meas(f_meas>0), S_xx_est(f_meas>0), 'r-','LineWidth',2,'DisplayName','S_{x,x} (model fit)')
     plot(f_meas(f_meas>0), S_zz_est(f_meas>0), 'c-','LineWidth',2,'DisplayName','S_{z,z} (model fit)')
     plot(f_meas(f_meas>0), S_xixi_est(f_meas>0), 'm-','LineWidth',2,'DisplayName','S_{\xi,\xi} (model fit)')
     
     box on
     set(gca,'XScale','log')
     set(gca,'YScale','log')
     xlim([xmin xmax])
     ylim([ymin ymax])
     title('DSH raw data')
     legend('Location','southwest');
     
   % smoothing
   %
     N_smooth = 1000;
     f_smooth    = 0;
     S_xx_smooth = 0;
     for i = 1 : N_smooth
       f_smooth    = f_smooth    + f(i:N_smooth:end);  
       S_xx_smooth = S_xx_smooth + S_xx(i:N_smooth:end);
     end
     f_smooth    = f_smooth/N_smooth;
     S_xx_smooth = S_xx_smooth/N_smooth;
     
     plot(f_smooth(f_smooth>0), S_xx_smooth(f_smooth>0), 'k-','LineWidth',3,'DisplayName','S_{x,x} (hidden, smoothed)')
    
    %}
     
     
     
     subplot(2,2,2); hold all     
     plot(f(f>0)          , S_xx(f>0), 'r-','LineWidth',1,'DisplayName','S_{x,x} (hidden)')
     plot(f_meas(f_meas>0), S_xhatxhat_inverse(f_meas>0), 'g-','LineWidth',1,'DisplayName','S_{x,x} (inv. reconst.)')
     plot(f_meas(f_meas>0), 1./Hsq(f_meas>0), 'k-','LineWidth',1,'DisplayName','|G|^2 (inv.)')
     plot(f_meas(f_meas>0), S_xixi_est(f_meas>0), 'm-','LineWidth',2,'DisplayName','S_{\xi,\xi} (model fit)')
     plot(f_meas(f_meas>0), S_xx_est(f_meas>0), 'r-','LineWidth',2,'DisplayName','S_{x,x} (model fit)')
     plot(f_meas(f_meas>0), SNR(f_meas>0), '-','Color',[1 0.7 0],'LineWidth',2,'DisplayName','SNR (model fit)')
     box on
     set(gca,'XScale','log')
     set(gca,'YScale','log')
     xlim([xmin xmax])
     ylim([ymin ymax])
     title('inverse filter')
     legend('Location','southwest');
     
     
     subplot(2,2,3); hold all
     plot(f(f>0)          , S_xx(f>0), 'r-','LineWidth',1,'DisplayName','S_{x,x} (hidden)')
     plot(f_meas(f_meas>0), S_xhatxhat_PSE(f_meas>0), 'g-','LineWidth',1,'DisplayName','S_{x,x} (PSE reconst.)')
     plot(f_meas(f_meas>0), G_PSE_sq(f_meas>0), 'k-','LineWidth',1,'DisplayName','|G|^2 (PSE)')
     plot(f_meas(f_meas>0), S_xixi_est(f_meas>0), 'm-','LineWidth',2,'DisplayName','S_{\xi,\xi} (model fit)')
     plot(f_meas(f_meas>0), S_xx_est(f_meas>0), 'r-','LineWidth',2,'DisplayName','S_{x,x} (model fit)')
     plot(f_meas(f_meas>0), SNR(f_meas>0), '-','Color',[1 0.7 0],'LineWidth',2,'DisplayName','SNR (model fit)')
     box on
     set(gca,'XScale','log')
     set(gca,'YScale','log')
     xlim([xmin xmax])
     ylim([ymin ymax])
     title('PSE filter')
     legend('Location','southwest');
     
     
     subplot(2,2,4); hold all
     plot(f(f>0)          , S_xx(f>0), 'r-','LineWidth',1,'DisplayName','S_{x,x} (hidden)')
     plot(f_meas(f_meas>0), S_xhatxhat_Wiener(f_meas>0), 'g-','LineWidth',1,'DisplayName','S_{x,x} (Wiener reconst.)')
     plot(f_meas(f_meas>0), G_Wiener_sq(f_meas>0), 'k-','LineWidth',1,'DisplayName','|G|^2 (Wiener)')
     plot(f_meas(f_meas>0), S_xixi_est(f_meas>0), 'm-','LineWidth',2,'DisplayName','S_{\xi,\xi} (model fit)')
     plot(f_meas(f_meas>0), S_xx_est(f_meas>0), 'r-','LineWidth',2,'DisplayName','S_{x,x} (model fit)')
     plot(f_meas(f_meas>0), SNR(f_meas>0), '-','Color',[1 0.7 0],'LineWidth',2,'DisplayName','SNR (model fit)')
     box on
     set(gca,'XScale','log')
     set(gca,'YScale','log') 
     xlim([xmin xmax])
     ylim([ymin ymax])
     title('Wiener filter')
     legend('Location','southwest');
     
     
     toc
     
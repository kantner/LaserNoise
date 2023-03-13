fprintf(1,'\n');
fprintf(1,'=== analytic PSDs of linearized system ===\n');

   % set pump currents for sweep
     I_range = [150 175 225 250 275 300 200]*1E-3; % the last one used in the next section for comparison with time-domain simulation          
     
   % set frequency points
     f_analytic = logspace(0,13,1000);     

   % analytic colored noise PSDs  S_FF = omega^-nu
     S_FF_P   = 1./((2*pi*f_analytic).^(par.nu(1)));
     S_FF_phi = 1./((2*pi*f_analytic).^(par.nu(2)));
     S_FF_N   = 1./((2*pi*f_analytic).^(par.nu(3)));     
   
   % allocate memory: white noise
     S_dPdP_WN         = zeros(size(f_analytic));
     S_dNdN_WN         = zeros(size(f_analytic));
     S_dphidphi_WN     = zeros(size(f_analytic));
     S_dPdN_WN         = zeros(size(f_analytic));
     S_dPdphi_WN       = zeros(size(f_analytic));
     S_dphidN_WN       = zeros(size(f_analytic));
     
   % allocate memory: colored noise  
     S_dPdP_CN_P       = zeros(size(f_analytic));
     S_dNdN_CN_P       = zeros(size(f_analytic));
     S_dphidphi_CN_P   = zeros(size(f_analytic));
     S_dPdN_CN_P       = zeros(size(f_analytic));
     S_dPdphi_CN_P     = zeros(size(f_analytic));
     S_dphidN_CN_P     = zeros(size(f_analytic));
     
     S_dPdP_CN_phi     = zeros(size(f_analytic));
     S_dNdN_CN_phi     = zeros(size(f_analytic));
     S_dphidphi_CN_phi = zeros(size(f_analytic));
     S_dPdN_CN_phi     = zeros(size(f_analytic));
     S_dPdphi_CN_phi   = zeros(size(f_analytic));
     S_dphidN_CN_phi   = zeros(size(f_analytic));
          
     S_dPdP_CN_N       = zeros(size(f_analytic));
     S_dNdN_CN_N       = zeros(size(f_analytic));
     S_dphidphi_CN_N   = zeros(size(f_analytic));
     S_dPdN_CN_N       = zeros(size(f_analytic));
     S_dPdphi_CN_N     = zeros(size(f_analytic));
     S_dphidN_CN_N     = zeros(size(f_analytic));
     
     P_stat = zeros(size(I_range));
     C_RIN  = zeros(size(I_range));
     K_RIN  = zeros(size(I_range));
     fc_RIN = zeros(size(I_range));
     C_FN   = zeros(size(I_range));
     K_FN   = zeros(size(I_range));
     fc_FN  = zeros(size(I_range));
   
     
   figure(2); clf; hold all;
   
   % loop over pump currents
   for i_pump = 1 : length(I_range)
       
     % set current 
       par.I = I_range(i_pump);
     
     % steady state values     
       [P_ss, N_ss] = steady_state(par);
    
     % compute system matrices  
       [J,G] = system_matrices_JG(par);
   
     % diffusion matrix (it holds 2*D = G*G' )
       D = 0.5 * G * (G'); 
     
     fprintf(1,'    compute PSDs at I = %.3f mA\n',par.I*1E3);
     tic   
   % loop over frequency   
     for i = 1 : length(f_analytic)    
     
     % invert system matrix  
       %K = - inverse(J + 1i*2*pi*f(i)*speye(3) );
       K = - (J + 1i*2*pi*f_analytic(i)*speye(3))\speye(3);
     
     % PSD matrices at selected frequency  
       S_WN     =   K * 2*D * (K');
       S_CN_P   =   K * sparse(1,1,1,3,3) * (K') * par.sigma(1)^2 *(2*P_ss)^2 * S_FF_P(i);
       S_CN_phi =   K * sparse(2,2,1,3,3) * (K') * par.sigma(2)^2 * S_FF_phi(i);     
       S_CN_N   =   K * sparse(3,3,1,3,3) * (K') * par.sigma(3)^2 * sqrt(N_ss) * S_FF_N(i);
     
     % pick out matrix components  
       S_dPdP_WN(i)     = S_WN(1,1);
       S_dphidphi_WN(i) = S_WN(2,2);
       S_dNdN_WN(i)     = S_WN(3,3);
       S_dPdN_WN(i)     = S_WN(1,3);
       S_dPdphi_WN(i)   = S_WN(1,2);
       S_dphidN_WN(i)   = S_WN(2,3);
     
       S_dPdP_CN_P(i)       = S_CN_P(1,1);
       S_dphidphi_CN_P(i)   = S_CN_P(2,2);
       S_dNdN_CN_P(i)       = S_CN_P(3,3);    
       S_dPdN_CN_P(i)       = S_CN_P(1,3);
       S_dPdphi_CN_P(i)     = S_CN_P(1,2);
       S_dphidN_CN_P(i)     = S_CN_P(2,3);
     
       S_dPdP_CN_phi(i)     = S_CN_phi(1,1);
       S_dphidphi_CN_phi(i) = S_CN_phi(2,2);
       S_dNdN_CN_phi(i)     = S_CN_phi(3,3);    
       S_dPdN_CN_phi(i)     = S_CN_phi(1,3);
       S_dPdphi_CN_phi(i)   = S_CN_phi(1,2);
       S_dphidN_CN_phi(i)   = S_CN_phi(2,3);     
          
       S_dPdP_CN_N(i)       = S_CN_N(1,1);
       S_dphidphi_CN_N(i)   = S_CN_N(2,2);
       S_dNdN_CN_N(i)       = S_CN_N(3,3);    
       S_dPdN_CN_N(i)       = S_CN_N(1,3);
       S_dPdphi_CN_N(i)     = S_CN_N(1,2);
       S_dphidN_CN_N(i)     = S_CN_N(2,3);     
     end
   
     S_dPdP     = real(S_dPdP_WN     + S_dPdP_CN_P     + S_dPdP_CN_phi     + S_dPdP_CN_N     );  
     S_dphidphi = real(S_dphidphi_WN + S_dphidphi_CN_P + S_dphidphi_CN_phi + S_dphidphi_CN_N );
     S_dNdN     = real(S_dNdN_WN     + S_dNdN_CN_P     + S_dNdN_CN_phi     + S_dNdN_CN_N     );
   
     S_dPdN     = S_dPdN_WN     + S_dPdN_CN_P     + S_dPdN_CN_phi     + S_dPdN_CN_N;
     S_dPdphi   = S_dPdphi_WN   + S_dPdphi_CN_P   + S_dPdphi_CN_phi   + S_dPdphi_CN_N;
     S_dphidN   = S_dphidN_WN   + S_dphidN_CN_P   + S_dphidN_CN_phi   + S_dphidN_CN_N;
   
   
     sgtitle('power spectral densities')
       
       subplot(2,2,1); hold all;     
         plot(f_analytic, S_dPdP_WN     * 1/(P_ss^2), 'c-.','LineWidth',2,'DisplayName','white noise')         
         plot(f_analytic, S_dPdP_CN_P   * 1/(P_ss^2), 'r:','LineWidth',2,'DisplayName','colored noise P')         
         plot(f_analytic, S_dPdP_CN_phi * 1/(P_ss^2), 'g:','LineWidth',2,'DisplayName','colored noise \phi')         
         plot(f_analytic, S_dPdP_CN_N   * 1/(P_ss^2), 'b:','LineWidth',2,'DisplayName','colored noise N')         
         plot(f_analytic, S_dPdP        * 1/(P_ss^2), 'k-','LineWidth',3,'DisplayName','white + colored noise')                  
         
         % least squares fit:  S_dPdP/(P_ss^2) = C/P^2 + (K/f)^nu         
           idx_range = find( (1E3 < f_analytic) & (f_analytic < 1E8) == 1);
           
           A = [1/P_ss^2,                                             mean(f_analytic(idx_range).^(-par.nu(2)));
                1/P_ss^2 * mean(f_analytic(idx_range).^(-par.nu(2))), mean(f_analytic(idx_range).^(-2*par.nu(2)))];
            
           b = [mean(S_dPdP(idx_range)/(P_ss^2)); mean(S_dPdP(idx_range)/(P_ss^2) .* 1./(f_analytic(idx_range).^par.nu(2)))];
           
           fit = A\b;
           
           P_stat(i_pump) = P_ss;
           C_RIN(i_pump)  = fit(1);
           K_RIN(i_pump)  = (fit(2))^(1/par.nu(2));           
           fc_RIN(i_pump) = K_RIN(i_pump) * (P_stat(i_pump)^2/C_RIN(i_pump))^(1/par.nu(2));
           
           plot(f_analytic(idx_range),     C_RIN(i_pump)/P_stat(i_pump)^2 + (K_RIN(i_pump)./f_analytic(idx_range)).^par.nu(2)   ,'m:','LineWidth',2,'DisplayName','least squares fit')
           plot(fc_RIN(i_pump),            C_RIN(i_pump)/P_stat(i_pump)^2 + (K_RIN(i_pump)./fc_RIN(i_pump)).^par.nu(2)          ,'ms','LineWidth',2,'DisplayName','corner frequency')
         
         box on
         axis tight
         set(gca,'XScale','log')
         set(gca,'YScale','log')
         xlabel('f [Hz]')
         ylabel('S_{\deltaP, \deltaP}/P_0^2 [Hz^2/Hz]')
         legend('Location','southwest','AutoUpdate','off')
         title('relative intensity noise')
         
         
       subplot(2,2,2); hold all;            
         plot(f_analytic, S_dphidphi_WN,     'c-.','LineWidth',2,'DisplayName','white noise')         
         plot(f_analytic, S_dphidphi_CN_P,   'r:','LineWidth',2,'DisplayName','colored noise P')         
         plot(f_analytic, S_dphidphi_CN_phi, 'g:','LineWidth',2,'DisplayName','colored noise \phi')         
         plot(f_analytic, S_dphidphi_CN_N,   'b:','LineWidth',2,'DisplayName','colored noise N')         
         plot(f_analytic, S_dphidphi,        'k-','LineWidth',3,'DisplayName','white + colored noise')         
         legend('Location','southwest','AutoUpdate','off')
         box on
         axis tight
         set(gca,'XScale','log')
         set(gca,'YScale','log')
         xlabel('f [Hz]')
         ylabel('S_{\delta\phi, \delta\phi} [Hz/Hz^2]')
         title('phase noise')
         
         
       subplot(2,2,3); hold all;     
         plot(f_analytic, S_dNdN_WN     * 1/(N_ss^2), 'c-.','LineWidth',2,'DisplayName','white noise')
         plot(f_analytic, S_dNdN_CN_P   * 1/(N_ss^2), 'r:','LineWidth',2,'DisplayName','colored noise P')
         plot(f_analytic, S_dNdN_CN_phi * 1/(N_ss^2), 'g:','LineWidth',2,'DisplayName','colored noise \phi')
         plot(f_analytic, S_dNdN_CN_N   * 1/(N_ss^2), 'b:','LineWidth',2,'DisplayName','colored noise N')
         plot(f_analytic, S_dNdN        * 1/(N_ss^2), 'k-','LineWidth',3,'DisplayName','white + colored noise')
         legend('Location','southwest','AutoUpdate','off')
         box on
         axis tight
         set(gca,'XScale','log')
         set(gca,'YScale','log')
         xlabel('f [Hz]')
         ylabel('S_{\deltaN, \deltaN}/N_0^2 [Hz^2/Hz]')
         title('relative carrier noise')

         
       subplot(2,2,4); hold all;            
         S_domegadomega        = (2*pi*f_analytic).^2 .* S_dphidphi;
         S_domegadomega_WN     = (2*pi*f_analytic).^2 .* S_dphidphi_WN;          
         S_domegadomega_CN_P   = (2*pi*f_analytic).^2 .* S_dphidphi_CN_P;
         S_domegadomega_CN_phi = (2*pi*f_analytic).^2 .* S_dphidphi_CN_phi;
         S_domegadomega_CN_N   = (2*pi*f_analytic).^2 .* S_dphidphi_CN_N;
         
       
         plot(f_analytic, S_domegadomega_WN,     'c-.','LineWidth',2,'DisplayName','white noise')
         plot(f_analytic, S_domegadomega_CN_P,   'r:','LineWidth',2,'DisplayName','colored noise P')
         plot(f_analytic, S_domegadomega_CN_phi, 'g:','LineWidth',2,'DisplayName','colored noise \phi')
         plot(f_analytic, S_domegadomega_CN_N,   'b:','LineWidth',2,'DisplayName','colored noise N')
         plot(f_analytic, S_domegadomega,        'k-','LineWidth',3,'DisplayName','white + colored noise')
         legend('Location','southwest','AutoUpdate','off')
         

         
        % least squares fit:  S_domegadomeg = C/P + (K/f)^nu
           idx_range = find( (1E3 < f_analytic) & (f_analytic < 1E8) == 1);
           
           A = [1/P_ss,                                             mean(f_analytic(idx_range).^(-par.nu(2)));
                1/P_ss * mean(f_analytic(idx_range).^(-par.nu(2))), mean(f_analytic(idx_range).^(-2*par.nu(2)))];
            
           b = [mean(S_domegadomega(idx_range)); mean(S_domegadomega(idx_range)./(f_analytic(idx_range).^par.nu(2)))];
           
           fit = A\b;
           
           C_FN(i_pump)   = fit(1);
           K_FN(i_pump)   = (fit(2))^(1/par.nu(2));
           P_stat(i_pump) = P_ss;
           fc_FN(i_pump)  = K_FN(i_pump) * (P_stat(i_pump)/C_FN(i_pump))^(1/par.nu(2));
           
           plot(f_analytic(idx_range),    C_FN(i_pump)/P_stat(i_pump) + (K_FN(i_pump)./f_analytic(idx_range)).^par.nu(2)   ,'m:','LineWidth',2,'DisplayName','least squares fit')
           plot(fc_FN(i_pump),            C_FN(i_pump)/P_stat(i_pump) + (K_FN(i_pump)./fc_FN(i_pump)).^par.nu(2)           ,'ms','LineWidth',2,'DisplayName','corner frequency')
         
           
         box on
         axis tight
         set(gca,'XScale','log')
         set(gca,'YScale','log')
         xlabel('f [Hz]')
         ylabel('S_{\delta\omega, \delta\omega} [Hz^2/Hz]')
         title('frequency noise')
      
         
         drawnow
         
    
         
    
   end
   
   
   
   
   
%% print fit parameters
   fprintf(1,'\n')
   fprintf(1,'    least squares fit of the RIN and FN PSD\n\n')
   fprintf(1,'    I [mA]\tP_stat\t\tC(RIN) [Hz]\tK(RIN) [Hz]\tfc(RIN) [Hz]\tC(FN) [Hz]\tK(FN) [Hz]\tfc(FN) [Hz]\n')
   for i_pump = 1 : length(I_range)
   fprintf(1,'    %.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n',I_range(i_pump)*1E3, P_stat(i_pump), C_RIN(i_pump), K_RIN(i_pump),fc_RIN(i_pump), C_FN(i_pump), K_FN(i_pump),fc_FN(i_pump));
   end
   fprintf(1,'\n');
   
   
   figure(3); clf; hold all;
     plot(P_stat, fc_RIN, 'ro-','LineWidth',2,'DisplayName','f_c (RIN)')
     plot(P_stat, fc_FN,  'go-','LineWidth',2,'DisplayName','f_c (FN)')
     box on
     xlabel('stationary photon number')
     ylabel('corner frequency [Hz]')
     axis tight
     title('corner frequency vs. stationary power')
     legend('location','northwest')
   
     
     
     
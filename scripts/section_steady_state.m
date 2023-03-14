   % compute threshold
     [N_th, I_th] = laserThreshold(par);

   % pump current range 
     I_range = logspace(log10(I_th*1E-2),log10(I_th*1E3),1001);
   
   % allocate memory
     N_bar   = zeros(length(I_range),1);
     P_bar   = zeros(length(I_range),1);
     
     g_bar        = zeros(length(I_range),1);     
     dg_dN_bar    = zeros(length(I_range),1);
     dg_dP_bar    = zeros(length(I_range),1);     
     g_sp_bar     = zeros(length(I_range),1);
     dg_sp_dN_bar = zeros(length(I_range),1);
     dg_sp_dP_bar = zeros(length(I_range),1);     
     r_bar        = zeros(length(I_range),1);          
     dr_dN_bar    = zeros(length(I_range),1);          
     
     Gamma_R_bar  = zeros(length(I_range),1);          
     Omega_R_bar  = zeros(length(I_range),1);          
     Gamma_P_bar  = zeros(length(I_range),1);          
     Gamma_N_bar  = zeros(length(I_range),1);          
     
     cov_PP = zeros(length(I_range),1);
     cov_PN = zeros(length(I_range),1);
     cov_NN = zeros(length(I_range),1);
     
   % sweep over pump current    
     for i = 1 : length(I_range)
       par.I = I_range(i);
     
     % compute steady stte
       [P_ss, N_ss] = steady_state(par);
       P_bar(i) = P_ss;
       N_bar(i) = N_ss;

     % steady state rates  
       [g_bar(i),    dg_dP_bar(i),    dg_dN_bar(i)]    = func_g(par, P_ss, N_ss);
       [g_sp_bar(i), dg_sp_dP_bar(i), dg_sp_dN_bar(i)] = func_g_sp(par, P_ss, N_ss);
       [r_bar(i),    dr_dN_bar(i)]                     = func_r(par, N_ss);     
     
     % relaxation rates
       [Gamma_R_bar(i), Omega_R_bar(i), Gamma_P_bar(i), Gamma_N_bar(i)] = relaxationRates(par);
     
     % steady state variance (for g(2) function)
       [J,G] = system_matrices_JG(par);
     
       % keep only cols and rows of the PN-subsystem (phase does not couple)
         J(2,:) = [];
         J(:,2) = [];
         G(2,:) = [];
         J = full(J);
         G = full(G);
       
       V = sylvester(J, J', -G*G');
       cov_PP(i) = V(1,1); % variance P
       cov_PN(i) = V(1,2); % variance N
       cov_NN(i) = V(2,2); % covariance P-N   
   end
  
   
 % plot steady state curves
   figure(1); clf; hold all;
   
   sgtitle('steady state characteristics and parameter dependencies')
   
   subplot(2,3,1); hold all;
     plot(I_range*1E3, P_bar, 'r-','LineWidth',2,'DisplayName','P')
     plot(I_range*1E3, N_bar, 'b-','LineWidth',2,'DisplayName','N')
     plot([1 1]*I_th*1E3, [min(P_bar),max(P_bar)],'r:','DisplayName','I_{th}')
     plot([min(I_range) max(I_range)]*1E3, [1 1]*N_th,'b:','DisplayName','N_{th}')     
     set(gca,'XScale','log')
     set(gca,'YScale','log')
     box on
     xlabel('pump current [mA]')
     ylabel('carrier number, photon number')
     title('steady state')
     axis tight
     legend('Location','southeast')

     
   subplot(2,3,2); hold all;     
     plot(I_range*1E3, par.Gamma*par.vg*g_bar,'k-.','LineWidth',2,'DisplayName','\Gamma v_g g')
     plot(I_range*1E3, par.Gamma*par.vg*g_bar.*P_bar,'k-','LineWidth',2,'DisplayName','\Gamma v_g g P')
     plot(I_range*1E3, par.Gamma*par.vg*g_sp_bar,'r-','LineWidth',2,'DisplayName','\Gamma v_g g_{sp}')
     plot(I_range*1E3, par.Gamma*par.vg*(g_sp_bar-g_bar),'b:','LineWidth',2,'DisplayName','\Gamma v_g g_{abs}')
     plot(I_range*1E3, r_bar,'g-','LineWidth',2,'DisplayName','r')
     plot(I_range*1E3, 1/par.tau_ph * ones(size(I_range)),'b-','LineWidth',2,'DisplayName','1/\tau_{ph}')
     plot(I_range*1E3, par.eta*I_range/elementaryCharge,'m--','LineWidth',2,'DisplayName','\eta I/q')
     plot(I_range*1E3, Gamma_P_bar,'r-.','LineWidth',2,'DisplayName','\Gamma_P')
     plot(I_range*1E3, Gamma_N_bar,'b-.','LineWidth',2,'DisplayName','\Gamma_N')
     plot(I_range*1E3, Gamma_R_bar,'k-.','LineWidth',4,'DisplayName','\Gamma_R')
     plot(I_range*1E3, Omega_R_bar,'c-.','LineWidth',4,'DisplayName','\Omega_R')
     box on
     xlabel('pump current [mA]')
     ylabel('rates [1/s]')
     set(gca,'XScale','log')
     set(gca,'YScale','log')
     legend('Location','northwest')
     title('rates and frequencies')
     axis tight
     
     
   N_range = linspace(0,5,101)*par.N_tr;
   
   subplot(2,3,3); hold all;     
     plot(N_range/par.N_tr, func_g(par,0,N_range),'k-','LineWidth',2,'DisplayName','g')
     plot(N_range/par.N_tr, func_g_sp(par,0,N_range),'r-','LineWidth',2,'DisplayName','g_{sp}')
     plot(N_range/par.N_tr, func_g_sp(par,0,N_range)-func_g(par,0,N_range),'b-','LineWidth',2,'DisplayName','g_{abs}')
     plot([min(N_range), max(N_range)]/par.N_tr, [0,0],'k--','LineWidth',1,'DisplayName','zero')
     box on
     xlabel('N/N_{tr}')
     ylabel('gain [1/m]')
     legend('Location','southeast')
     title('gain coefficient')
     axis tight
     
   
   subplot(2,3,4); hold all;
     g2      = @(VarN,N) (VarN + N.*(N-1))./(N.*N);
     g2_eval = g2(cov_PP,P_bar);
     plot(I_range*1E3, g2_eval,'r-','LineWidth',2,'DisplayName','g^{(2)}(0)')
     box on
     xlabel('pump current [mA]')
     ylabel('g^{(2)}(0)')   
     set(gca,'XScale','log')
     axis tight
     title('second-order correlation funtion')     
     plot([1 1]*I_th*1E3, [0,2.5],'r:','DisplayName','threshold')
     plot(I_range*1E3, 1*ones(size(I_range)),'k--','LineWidth',2,'DisplayName','coherent')
     plot(I_range*1E3, 2*ones(size(I_range)),'k:','LineWidth',2,'DisplayName','thermal')
     legend('Location','northeast')
     ylim([0.5,2.5])
     
     
   subplot(2,3,5); hold all;
     plot(I_range*1E3, cov_PP./P_bar,'r-','LineWidth',2,'DisplayName', 'Var(P)/P (coherent)')
     plot(I_range*1E3, cov_PP./(P_bar.*(P_bar+1)),'k--','LineWidth',2,'DisplayName', 'Var(P)/(P(P+1)) (thermal)')
     plot(I_range*1E3, cov_NN./N_bar,'b-','LineWidth',2,'DisplayName', 'Var(N)/N')
     box on
     xlabel('I [mA]')
     ylabel('Var/Mean')   
     set(gca,'XScale','log')
     set(gca,'YScale','log')
     axis tight
     title('variance')          
     plot([1 1]*I_th*1E3, [min(cov_PP./(P_bar.*(P_bar+1))),max(cov_PP./(P_bar.*(P_bar+1)))],'r:','DisplayName', 'threshold')
     legend('Location','southwest')
     
  
   subplot(2,3,6); hold all;     
     plot(N_range/par.N_tr, func_beta_sp(par,0,N_range),'g-','LineWidth',2,'DisplayName', '\beta_{sp}(N)')
     box on
     xlabel('N/N_{tr}')
     ylabel('\beta_{sp}')
     axis tight
     title('spontaneous emission \beta_{sp}-factor')  
     legend('Location','northeast')



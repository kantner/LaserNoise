function [x] = correct_phase_jumps(x)

  % Routine to correct a function that exhibits articifical jumps +/- pi
  % stemming from branch cuts in an inverse trigonometric function.
  % Increments are expected to be actually small, such that jumps close to
  % +/- pi are regarded as artifacts that are corrected in this routine.
  %
  % INPUT
  %   x
  % OUTPUT
  %   x_corr
  
  DEBUG = 0; % switch on plotting
  
    
  %% step 1: histogram of increments
   % search for increments   
     
     
   
     [vals, edges] = histcounts(diff(x),'Normalization','pdf');
     
     
     if DEBUG == 1
     figure(1001); clf; hold all;              
       h = histogram(diff(x),'Normalization','pdf');
       box on;
       title('histogram of increments')       
       xlabel('increment')
       ylabel('counts (normalized)')
     end
      
      
      
    
    % compute mean and variance
      x_axis     = 0.5*(edges(2:end) +  edges(1:end-1));
      mean_value = trapz(x_axis, x_axis.*vals);
      variance   = trapz(x_axis, x_axis.^2 .* vals) - mean_value^2;
    
    % set threshold  
      threshold = 5 * sqrt(variance);
      
      if DEBUG == 1
      gaussian_pdf = @(x, mu, sigma_sq) exp(-0.5*(x-mu).^2/sigma_sq)/sqrt(2*pi*sigma_sq);
      
      plot(x_axis, gaussian_pdf(x_axis, mean_value, variance), 'r-')
      
      plot( [1 1]*(+1)*threshold, [0 max(h.Values)], 'r--')
      plot( [1 1]*(-1)*threshold, [0 max(h.Values)], 'r--')
      end
      
      

      
     
  %% step 2: correction  --- slow method
  %{
   % identify points where index jumps occur
     idx_up   = find(x(1:end-1)-x(2:end) > +threshold);
     idx_down = find(x(1:end-1)-x(2:end) < -threshold);
     
     fprintf(1,'  identified %d with jumps close to +pi\n',length(idx_up))
     fprintf(1,'  identified %d with jumps close to -pi\n',length(idx_down))
     
     if DEBUG == 1
     figure(1002);clf; hold all;
     plot(1:length(x),x          ,'k-','DisplayName','x')
     plot(idx_up,     x(idx_up)  ,'ro','DisplayName','phase jump +\pi')
     plot(idx_down,   x(idx_down),'go','DisplayName','phase jump -\pi')
     box on
     xlabel('step index')
     ylabel('x')
     legend()
     end
     
   % correction 
     tic
     for i = 1 : length(idx_up)
       x(idx_up(i)+1:end)   = x(idx_up(i)+1:end)   + pi;       
     end     
     for i = 1 : length(idx_down)
       x(idx_down(i)+1:end) = x(idx_down(i)+1:end) - pi;
     end
     toc
     
     if DEBUG == 1
     plot(1:length(x),x,'b-','DisplayName','x (corrected - slow)')
     end
  %}
     
  %% step 2: correction --- fast method
    
     if DEBUG == 1
       figure(1002);clf; hold all;
         plot(1:length(x),x          ,'k-','DisplayName','x')
         box on
         xlabel('step index')
         ylabel('x')
         legend()     
     end
  
     %tic
     add_pi      = cumsum(x(1:end-1)-x(2:end) > +threshold)*pi;
     subtract_pi = cumsum(x(1:end-1)-x(2:end) < -threshold)*pi;
     
     x(2:end)    = x(2:end) + add_pi - subtract_pi;
     %toc
     
     if DEBUG == 1
         plot(1:length(x),x          ,'m--','DisplayName','x (corrected - fast)')     
     end
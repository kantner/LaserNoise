function [x,success] = newton(f,df,x0,maxIterations,tol,contr_threshold,iroundmax,damp_mode,quiet_mode,newton_mode)



% NEWTON ITERATION
% input: f(x), df(x), x0 (initial guess)

% default
%{
  maxIterations   = 100;
  tol             = 1E-36;
  contr_threshold = 1E-5;
  iroundmax       = 5;  
%}  
  %newton_mode = 0;
  % 0: vectorial (standard: solve a system of equations with a single guess)
  % 1: element-wise (solve a single equation with multiple initial guesses in parallel)

  %damp_mode   = 0;
  % 0: no damping
  % 1: use residuum

  %quiet_mode = 0;
  % 0: print results
  % 1: silence = print no results


% user-specified: overwrite default values
%{
  if nargin >= 4
    maxIterations   = varargin{1};
  end
  if nargin >= 5
    tol             = varargin{2};
  end
  if nargin >= 6
    contr_threshold = varargin{3};
  end
  if nargin >= 7
    iroundmax       = varargin{4};
  end    
  if nargin >= 8
    newton_mode     = varargin{5};
  end    
  if nargin >= 9
    damp_mode       = varargin{6};
  end   
 %}
  
% control 
  continue_iteration = 1;
  iter               = 0;
  count_roundoff     = 0; 
  iround             = 0;



% norm
  norm_max = @(x) max(max(abs(x))); % default
  norm_L2  = @(x) sqrt(x'*x);

  norm = @(x) norm_max(x);



x = x0;

while continue_iteration == 1
   
  % iteration counter
    iter = iter + 1;
    
  % compute update
    switch(newton_mode)
      case 0 % vectorial
        if iter == 1 && quiet_mode == 0            
        fprintf(1,'Newton method: vectorial mode\n')
        end    
        dx  = - df(x)\f(x);
      case 1 % element-wise
        if iter == 1 && quiet_mode == 0
        fprintf(1,'Newton method: element-wise mode\n')
        end    
        dx  = - df(x).\f(x);
    end

      
    if iter == 1 && quiet_mode == 0
    fprintf(1,'i\tres\t\tcontr\t\tdamp\t\tiround\n')
    end
    
        
  % compute residuum  
    if iter > 1
    res_old = res;    
    end    
    res = norm(dx);
  
    
  % compute damping factor
    switch(damp_mode)
      case 0 % no damping
        damp = 1;
      case 1 % use residuum
         damp_min = 0.1; 
         damp = max( 1/(1+res), damp_min);
    end
    
  % update
    x = x + damp * dx;
    
    
  % contraction
    if iter > 1
    contr = res/res_old;
    end    
  
  % print  
    if quiet_mode == 0
    if iter == 1
      fprintf(1,'%d\t%.4e\t\t\t%.4e\n',iter, res, damp);
    elseif count_roundoff == 1
      fprintf(1,'%d\t%.4e\t%.4e\t%.4e\t%d\n',iter, res, contr, damp, iround);
    else
      fprintf(1,'%d\t%.4e\t%.4e\t%.4e\n',iter, res, contr, damp);
    end
    end
    
    
    
    %{
    figure(1);hold all;
    plot([x-dx x],[f(x-dx) f(x)],'bo-')
    drawnow
    %}
    
    
  % control
    if iter >= maxIterations
      if quiet_mode == 0
      fprintf('Newton iteration failed: max number of iterations reached\n')
      end
      success = 0;
      continue_iteration = 0;
    end
    
    if res < tol
      if quiet_mode == 0
      fprintf('Newton iteration successful: Tolerance reached\n')
      end
      success = 1;
      continue_iteration = 0; 
    end
    
    if iter > 1
    if contr < contr_threshold && count_roundoff == 0 && iter > 1
      % we should at least have completed 2 iterations before claiming this <<< --- no! set to "iter > 0"
      count_roundoff = 1;      
    end
    end
    
    if count_roundoff == 1
      iround = iround + 1; 
    end
    
    if iround >= iroundmax
      if quiet_mode == 0
      fprintf('Newton iteration successful: Round-off error reached.\n')
      end
      success = 1;
      continue_iteration = 0;     
    end
    
    %{
    if continue_iteration == 0
    if success == 1
       figure(1);hold all;
        plot(x,fx,'gs') 
    end
    end
    %}
    
    
    
    
    

    
end


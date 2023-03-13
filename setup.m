%% parameters
   physical_constants
   

 % device geometry
   par.d       = 5*10E-9; % active region thickness | unit m
   par.W       =    5E-6; % active region width     | unit m
   par.L       =    5E-3; % active region length    | unit m
 
 % wavelength
   par.lambda0     = 1064E-9; % vacuum wavelength     | unit m
   
 % temperature  
   par.temperature = 300;     % temperature  | unit K
  
 % light-matter interaction
   par.ng      = 3.9;       % group index                  | unit 1
   par.g0      = 354000;    % gain coefficient             | unit 1/m   
   par.n_tr    = 2E18*1E6;  % transparency carrier density | unit m^-3      
   par.Gamma   = 0.01;      % opt. confinement factor      | unit 1
   par.eps     = 1E-8;      % gain compression coefficient | unit 1
   par.alpha_H = 3.0;       % Henry factor                 | unit 1
      
 % recombination
   par.A       = 1*1E8;    % SRH recombination rate         | unit s^-1
   par.B       = 1*1E-16;  % bimolecular recombination rate | unit m^3 s^-1
   par.C       = 1*4E-42;  % Auger recombination rate       | unit m^6 s^-1

 % injection
   par.I       = 200E-3; % injection current    | unit A
   par.eta     = 0.9;    % injection efficiency | unit 1 
   
 % photon lifetime
   par.tau_ph = 5E-10;

 % dependent parameters
   par.vg          = vacuumSpeedOfLight/par.ng;                    % group velocity               | unit m/s
   par.V           = par.d * par.W * par.L;                        % active region volume         | unit m^3
   par.N_tr        = par.V * par.n_tr;                             % transparency carrier number  | unit 1
   par.omega0      = 2*pi*vacuumSpeedOfLight/par.lambda0;          % central (angular) frequency  | unit 1/s
   par.P_th        = 1/(exp(reducedPlanckConstant*par.omega0/(boltzmannConstant*par.temperature)) - 1); % thermal equilibrium photon number | unit 1
   
 % colored noise parameters for (F_P F_phi F_N) 
   par.nu      = [1.4, 1.4, 1.0];          % power law coefficients
   par.N_F     = [1,   1,   1]    *20;     % number of fluctuators for discretization
   par.sigma   = [5E5, 5E5, 1E9];          % amplitudes
   
 % SDH measurement
   par.tau_d      = 10E-6;  % delay time                 | unit s
   par.eta_det    = 1;      % detector efficiency        | unit ????   
   par.deltaOmega = 100E6;  % frequency offset           | unit 1/s
   par.sigma_IQ   = 2E3;    % I-Q  white noise level     | unit 1/s
   par.r0         = 0.7;
   par.eta_out    = (1-par.r0^2)*par.vg*par.omega0*reducedPlanckConstant/par.L;   % conversion factor to power | unit W
   
 % Hann windowing
   par.hann_window_switch = 1;
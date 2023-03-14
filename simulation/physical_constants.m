%% some physical constants
% all in SI basis unit

  planckConstant        = 6.62607015 * 1E-34;    % (h) exact                       
  reducedPlanckConstant = planckConstant/(2.0*pi); % (hbar) exact
  elementaryCharge      = 1.602176634 * 1E-19;   % exact
  boltzmannConstant     = 1.380649 * 1E-23;      % exact                        
  vacuumSpeedOfLight    = 299792458;             % exact
  
  vacuumDielectricConstant = 8.8541878128 * 1E-12;  
  magneticConstant         = 1.0/(vacuumSpeedOfLight*vacuumSpeedOfLight*vacuumDielectricConstant);
                             
  electronMass       = 9.1093837015 * 1E-31;
                       
  bohrRadius         = 4.0*pi*vacuumDielectricConstant*reducedPlanckConstant*reducedPlanckConstant/(electronMass*elementaryCharge*elementaryCharge);
  bohrMagneton       = elementaryCharge*reducedPlanckConstant/(2.0*electronMass);

  
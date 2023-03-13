clear all
clc
 
 % add subdirectories to path
   addpath('distributions/')
   addpath('functions/')
      
 % suppress warnings
   warning('off','MATLAB:Axes:NegativeDataInLogAxis')
   warning('off','MATLAB:plot:IgnoreImaginaryXYPart')
   
   
%% initialize random number generator
   rng(18199,'twister')
   
%% set parameters
   setup;
   
%% steady states
   section_steady_state;
  
%% noise power spectral densities
   section_PSDs;
         
%% time domain simulation
   section_time_domain_simulation;
            
%% self-heterodyne detection
   section_SDH_measurement;

%% Wiener filtering of detected signals   
   section_filtering
   
   
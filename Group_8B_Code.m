% MECE 301 - Engineering Applications Lab
% Final Project - Pneumatic Tube
% Written by Jacob, Jason, Jonah, and Tim

% Clear command window, variables, and figures
clc; clear; close all;

%% PHYSICAL CONSTANTS
% USE FORMAT : [estimate_value uncertainty]

% Room Temp and Atmospheric Pressure Air
T_room = [19 3]; % C
p_atm = [101.3 2]; % kPa

% Dimensions
pipe_diam = [10 1]; % m
carrier_diam = [10 1]; % m

carrier_mass = [10 1]; % kg

g = 9.81; % m/s^2


%% FUNCTIONS
% Calculate density of air as a function of temperature and pressure using ideal gas law
% INPUTS:
%  * pres_Pa = absolute pressure in Pascals (N/m^2)
%  * temp_C  = temperature in Celcius
% OUTPUTS:
%  * density = density of dry air in kg/m^3
function [density] = calcAirDensity(pres_Pa, temp_C)
    % p*V = m*R_specific*T --> m/V = p / (R_specific * T)
    R_specific = 287.05; % J / (kg * K) - for dry air
    density = pres_Pa ./ (R_specific .* (temp_C + 273.15));
end

% INPUTS :
%  * length : length of the carrier [m]
%  * sys_params : all relevant physical parameters of the launcher
%       1 - 
%  * printing : false = skip plotting, true = plot graph of time dependent terms
% OUTPUT : terminal velocity of ball (in m/s^2)
function [term_vel] = modelLaunch(length, sys_params, printing)

% GIVEN CONSTANTS - taken from input array
    
% CALCULATED CONSTANTS  
    
% TIME DEPENDENT VARIABLES
    n_steps = 500; % Number of steps in the loop
    delta_t = 0.01;  % Amount of time between each step

    time = zeros(n_steps,1);
    position = zeros(n_steps,1);
    velocity = zeros(n_steps,1);
    accel = zeros(n_steps,1);
    
% INITIAL CONDITIONS
     
% MAIN EULER CALCULATION LOOP
    n = 1;
    while n < n_steps

        % Increment step counter
        n = n + 1;
    end

% SAMPLE PLOT PRINTING
    if printing
        
    end
end
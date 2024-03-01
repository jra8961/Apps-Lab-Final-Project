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

%% RUN SIMULATION
% Number of points on the length vs. velocity graph
n_pts = 20;

% Setup arrays for storing carrier length and terminal velocity
terminal_velocity = zeros(n_pts,3);
carrier_length = linspace(1,10,n_pts); % kPa

% Store all physical parameters in an array to make it easy to pass into function
system_parameters = [ ];

% Store all variations in system parameters for easy access during loop
system_variations = [ ];

% Loop for each carrier length
for m = 1:n_pts
    % Calculate estimate exit velocity for that vacuum pressure
    terminal_velocity(m,2) = modelTube(carrier_length,system_parameters,false);
    
    % Initialize variables for storing variations caused by each parameter
    max_error = 0;
    min_error = 0;
    
    % For each system parameter ...
    for n = 1:length(system_parameters)
        % Create temporary copy of parameters
        temp_params = system_parameters;
        
        % Calculate the MIN variation for one parameter and calculate the effect it has
        temp_params(n) = system_parameters(n) + system_variations(n);
        temp_error_1 = modelLaunch(vacuum(m),temp_params,false) - terminal_velocity(m,2);

        % Calculate the MIN variation for one parameter and calculate the effect it has
        temp_params(n) = system_parameters(n) - system_variations(n);
        temp_error_2 = modelLaunch(vacuum(m),temp_params,false) - terminal_velocity(m,2);
        
        % One variation will cause a positive error and the other will create a negative error
        % Determine which is positive and add to corrisponding error sum
        % Values are squared here for RSS combination of errors. Square root is later in code
        if (temp_error_1 > 0 && temp_error_2 < 0)
            max_error = max_error + temp_error_1^2;
            min_error = max_error + temp_error_2^2;
        elseif (temp_error_1 < 0 && temp_error_2 > 0)
            max_error = max_error + temp_error_2^2;
            min_error = max_error + temp_error_1^2;
        end
    end

    % Calculate the error bounds for this vacuum pressure
    % Apply square root of RSS combination here
    terminal_velocity(m,1) = terminal_velocity(m,2) - sqrt(min_error);
    terminal_velocity(m,3) = terminal_velocity(m,2) + sqrt(max_error);
end

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
%  * printing : false = skip plotting, true = plot sample graph of time dependent terms
% OUTPUT : terminal velocity of ball (in m/s^2)
function [term_vel] = modelTube(length, sys_params, printing)

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

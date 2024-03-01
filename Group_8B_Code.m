% MECE 301 - Engineering Applications Lab
% Final Project - Pneumatic Tube
% Written by Jacob, Jason, Jonah, and Tim

% Clear command window, variables, and figures
clc; clear; close all;

%% PHYSICAL CONSTANTS
% USE FORMAT : [estimate_value uncertainty]

% Room Temp and Atmospheric Pressure Air
room_temp = [19 3]; % C
p_atm = [101.3 2] * 10^3; % Pa
room_density = calcAirDensity(p_atm, room_temp);

% ALL PARAMETERS BELOW ARE PLACE HOLDERS !!!!!
% PLEASE DO RESEARCH AND FIND ACTUAL VALUES !!!!!

% Dynamic viscosity of air
viscosity = 100;

% Tube Dimensions
tube_diam = [10 1]; % m
tube_length = [60 1]; % m

% Carrier Dimensions
carrier_diam = [10 1]; % m
carrier_mass = [10 1]; % kg
min_carrier_length = 10; % m
max_carrier_length = 10; % m

% Gap between the carrier and tube
air_gap = (tube_diam - pipe_diam) / 2;

% Gravitational acceleration
g = 9.81; % m/s^2

%% RUN SIMULATION
% Number of points on the length vs. velocity graph
n_pts = 20;

% Setup arrays for storing carrier length and terminal velocity
terminal_velocity = zeros(n_pts,3);
carrier_length = linspace(min_carrier_length,max_carrier_length,n_pts); % kPa

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
        
        % Calculate the MAX variation for one parameter and calculate the effect it has
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

% Use Euler's Method to solve for velocity of the carrier over time
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
    carrier_velocity = zeros(n_steps,1);
    avg_air_velocity = zeros(n_steps,1);
    accel = zeros(n_steps,1);

    p_inside = zeros(n_steps,1);
    inside_volume = zeros(n_steps,1);
    inside_mass = zeros(n_steps,1);
    inside_density = zeros(n_steps,1);
    
% INITIAL CONDITIONS (any variables that are not expicitly set here are initially zero)
    p_inside(1) = p_atm;
    inside_volume(1) = tube_length * pi() * tube_diam^2 /4;
    inside_mass(1) = inside_volume(1) * room_density;
    inside_density = inside_volume(1) / inside_mass(1);

% MAIN EULER CALCULATION LOOP
    n = 1;
    while n < n_steps
        % Calculate new current time
        time(n+1) = delta_t * n;

        % Determine if the entire carrier is inside the tube
        % If not, only count the length of the carrier below the top of the tube
        inside_length = length;
        if(position(n) < length)
            inside_length = position(n)
        end

        % Pressure force between air above and below carrier
        f_pressure = (p_atm-p_inside) * pi()*tube_diam^2/4;

        % Shear force along sides of carrier
        f_shear = 0;

        % Resulting acceleration from gravity, pressure, and shear
        accel(n) = g + (f_pressure - f_shear) / carrier_mass;

        % Resulting velocity and position
        carrier_velocity(n+1) = carrier_velocity(n) + accel(n) * delta_t;
        position(n+1) = position(n) + carrier_velocity(n) * delta_t;
        
        % Resulting average velocity in annulus
% !!!!! DOUBLE CHECK THIS DERIVATION !!!!!!
        pressure_gradient = (p_atm-p_inside) / inside_length;
        avg_air_velocity(n+1) = pressure_gradient / (2*viscosity) * air_gap^2 * (air_gap/3 - 0.5) + ...
            carrier_velocity*air_gap/2;

        % Air displaced by the carrier moving down
        vol_displaced = (position(n+1) - position(n)) * pi()*carrier_diam^2/4;
        mass_displaced = vol_displaced * inside_density(n);
        
        % How much mass left through the annulus?
        % mass flow rate = density * velocity * area
        exited_mass = inside_density(n) * avg_air_velocity * pi()/4 (tube_diam^2 - carrier_diam^2) * delta_t
        
        % How much air is left underneath the carrier?
        inside_mass(n+1) = inside_mass(n) - exited_mass;
        inside_volume(n+1) = (position(n+1) - position(n)) * pi()*tube_diam^2/4;
        inside_density(n+1) = inside_mass(n+1) / inside_mass(n+1);

        % Use ideal gas law to determine what the next pressure
        p_inside(n+1) = inside_volume(1)/inside_volume(n+1) * inside_mass(1)/inside_mass(n+1) * p_inside(1);

        
        

        % Increment step counter
        n = n + 1;
    end

% SAMPLE PLOT PRINTING
    if printing
        
    end
end
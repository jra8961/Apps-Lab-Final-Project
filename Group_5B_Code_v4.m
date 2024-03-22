% MECE 301 - Engineering Applications Lab
% Final Project - Pneumatic Tube
% Written by Jacob, Jason, Jonah, and Tim

% Clear command window, variables, and figures
clc; clear; close all;

%% PHYSICAL CONSTANTS
% USE FORMAT : [estimate_value uncertainty]

% Sea-level density and Atmospheric Pressure
%https://www.earthdata.nasa.gov/topics/atmosphere/atmospheric-pressure/air-mass-density#:~:text=Pure%2C%20dry%20air%20has%20a,a%20pressure%20of%20101.325%20kPa.
air_density = [1.293 0.01]; % kg/m3
p_atmos = [101.3 1] * 10^3; % Pa

% Dynamic viscosity of air
% SOURCE: https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html
air_viscosity = [1.822 0.0001822] * 10^-5; %N*s/m2

% Tube Dimensions
tube_diameter = [1.772 0.003] * 0.0254; % in -> m
tube_length = [48 0.5] * 0.0254; % in -> m

% Carrier Dimensions
carrier_diameter = [1.735 0.002] * 0.0254; % in -> m
carrier_mass = [0.250 0.005]; % kg
min_carrier_length = 0.04; % m
max_carrier_length = 0.16; % m

%% RUN SIMULATION
% Number of points on the length vs. velocity graph
n_pts = 16;

% Setup arrays for storing carrier length and terminal velocity
terminal_velocity = zeros(n_pts,3);
carrier_length = linspace(min_carrier_length,max_carrier_length,n_pts); % kPa

% Store all physical parameters in an array to make it easy to pass into function
system_parameters = [air_density(1), p_atmos(1), air_viscosity(1), tube_diameter(1), tube_length(1), ...
    carrier_diameter(1), carrier_mass(1)];

% Store all variations in system parameters for easy access during loop
system_variations = [air_density(2), p_atmos(2), air_viscosity(2), tube_diameter(2), tube_length(2), ...
    carrier_diameter(2), carrier_mass(2)];

% Loop for each carrier length
for m = 1:n_pts
    % Calculate estimate exit velocity for that vacuum pressure
    terminal_velocity(m,2) = modelTube(carrier_length(m),system_parameters,false);
    
    % Initialize variables for storing variations caused by each parameter
    upper_bound_error = 0;
    lower_bound_error = 0;
    
    % For each system parameter ...
    for n = 1:length(system_parameters)
        % Create temporary copy of parameters
        temp_params = system_parameters;
        
        % Calculate the MAX variation for one parameter and calculate the effect it has
        temp_params(n) = system_parameters(n) + system_variations(n);
        temp_error_1 = modelTube(carrier_length(m),temp_params,false) - terminal_velocity(m,2);

        % Calculate the MIN variation for one parameter and calculate the effect it has
        temp_params(n) = system_parameters(n) - system_variations(n);
        temp_error_2 = modelTube(carrier_length(m),temp_params,false) - terminal_velocity(m,2);
        
        % Check which variation make the velocity estimate bigger and which
        % one makes it smaller (-error = larger velocity bc y-axis points upward)
        if (temp_error_1) > 0 && (temp_error_2) < 0
            upper_bound_error = upper_bound_error + temp_error_1^2;
            lower_bound_error = lower_bound_error + temp_error_2^2;
        else
            upper_bound_error = upper_bound_error + temp_error_2^2;
            lower_bound_error = lower_bound_error + temp_error_1^2;
        end
    end

    % Calculate the error bounds for this vacuum pressure
    % Apply square root of RSS combination here
    terminal_velocity(m,1) = terminal_velocity(m,2) - sqrt(lower_bound_error);
    terminal_velocity(m,3) = terminal_velocity(m,2) + sqrt(upper_bound_error);
end

%% PLOTTING
% Sample plot of carrier position, velocity, and acceleration
modelTube(max_carrier_length,system_parameters,true);

% Plot of carrier length vs. terminal velocity
figure();
hold on
% Estimate velocity
plot(carrier_length,terminal_velocity(:,2));
% Error bounds of estimate
plot(carrier_length,terminal_velocity(:,1),'--');
plot(carrier_length,terminal_velocity(:,3),'--');
hold off

% Make figure look pretty
grid on;
xlabel("Carrier Length [m]"); ylabel("Terminal Velocity [m/s]");
legend("Estimated Results", "Min Uncertainty Bound", "Max Uncertainty Bound");
%ylim([0 0.2]);
legend("Location","southeast");

%% FUNCTIONS
% Use Euler's Method to solve for velocity of the carrier over time
% INPUTS :
%  * length : length of the carrier [m]
%  * sys_params : all relevant physical parameters of the launcher
%       1 - room_density : density of dry, sea-level air [kg/m^3]
%       2 - p_atm : atmospheric pressure [in Pa or N/m^2]
%       3 - viscosity : dynamic visocity of air [in N*s/m^2]
%       4 - tube_diam : inner diameter of the tube [in m]
%       5 - tube_len : length of the tube [in m]
%       6 - carrier_diam : outer diameter of the carrier [in m]
%       7 - carrier_mass : mass of the carrier [in kg]
%  * printing : false = skip plotting, true = plot sample graph of time dependent terms
% OUTPUT : 
%  * term_vel : terminal velocity of ball [in m/s]
function [term_vel] = modelTube(carrier_len, sys_params, printing)

% GIVEN CONSTANTS
   % Gravitational acceleration
    g = 9.81; % m/s^2 
    
    % Extract parameters from input vector
    room_density = sys_params(1);
    p_atm = sys_params(2);
    viscosity = sys_params(3);
    tube_diam = sys_params(4);
    tube_len = sys_params(5);
    carrier_diam = sys_params(6);
    carrier_mass = sys_params(7);

% CALCULATED CONSTANTS  
    % Gap between the carrier and tube
    air_gap = (tube_diam - carrier_diam) / 2;

% TIME DEPENDENT VARIABLES
    n_steps = 5000; % Number of steps in the loop
    delta_t = 0.00006;  % Amount of time between each step

    time = zeros(n_steps,1);
    position = zeros(n_steps,1);
    carrier_velocity = zeros(n_steps,1);
    avg_air_velocity = zeros(n_steps,1);
    accel = zeros(n_steps,1);
    
    pressure_gradient = zeros(n_steps,1);

    p_inside = zeros(n_steps,1);
    inside_volume = zeros(n_steps,1);
    inside_mass = zeros(n_steps,1);
    inside_density = zeros(n_steps,1);
    
% INITIAL CONDITIONS (any variables that are not expicitly set here are initially zero)
    accel(1) = -g;
    p_inside(1) = p_atm;
    inside_volume(1) = tube_len * pi()/4 * tube_diam^2;
    inside_mass(1) = inside_volume(1) * room_density;
    inside_density(1) = room_density;

% Percent change in velocity loop-to-loop for determining stability
    percent_change = 1;


% MAIN EULER CALCULATION LOOP
    n = 1;
    while (percent_change > 10^-11) && (abs(position(n)) < tube_len)
        % Calculate new current time
        time(n+1) = delta_t * n;

        % Pressure force between air above and below carrier
        f_pressure = (p_inside(n) - p_atm) * pi()/4*carrier_diam^2;

        % Pressure gradient in the annulus
        pressure_gradient(n) = (p_atm - p_inside(n)) / carrier_len;
        
        % Derivative of the air velocity profile at the face of the carrier    
% !!!!! DOUBLE CHECK THIS DERIVATION !!!!!!
        air_velocity_derivative = -carrier_velocity(n)/air_gap - air_gap*pressure_gradient(n)/(2*viscosity);

        % Shear force acting on the side face of the carrier
        f_shear = viscosity * air_velocity_derivative * carrier_len*pi()*carrier_diam;

        % Resulting acceleration from gravity, pressure, and shear
        accel(n+1) = -g + (f_pressure + f_shear) / carrier_mass;

        % Resulting velocity and position
        carrier_velocity(n+1) = carrier_velocity(n) + accel(n) * delta_t;
        position(n+1) = position(n) + carrier_velocity(n) * delta_t;
        
        % Resulting average velocity in annulus
% !!!!! DOUBLE CHECK THIS DERIVATION !!!!!!
        avg_air_velocity(n+1) = carrier_velocity(n)/2 - pressure_gradient(n)/(12*viscosity)*air_gap^2;

        % Air displaced by the carrier moving down
        vol_displaced = abs(position(n+1) - position(n)) * pi()/4*carrier_diam^2;
        mass_displaced = vol_displaced * inside_density(n);
        
        % How much mass left through the annulus?
        % mass flow rate = density * velocity * area
        exited_mass = inside_density(n) * avg_air_velocity(n) * pi()/4 * (tube_diam^2 - carrier_diam^2) * delta_t;

        % How much air is left underneath the carrier?
        inside_mass(n+1) = inside_mass(n) - exited_mass;
        inside_volume(n+1) = (tube_len + position(n)) * pi()/4*tube_diam^2;
        inside_density(n+1) = inside_mass(n+1) / inside_volume(n+1);

        % Use ideal gas law to determine what the next pressure
        p_inside(n+1) = inside_volume(1)/inside_volume(n) * inside_mass(n)/inside_mass(1) * p_inside(1);

        % Determine how much velocity has changed from last loop
        if (n > 1 && carrier_velocity(n-1) ~= 0)
            percent_change = abs( (carrier_velocity(n+1) - carrier_velocity(n))/carrier_velocity(n) );
        end

        % Increment step counter
        n = n + 1;
    end
    
    % The last velocity value (where the velocity settled) is the terminal velocity
    term_vel = carrier_velocity(end);

% SAMPLE PLOT PRINTING
    if printing
        figure();
        
        % Plot velocity of the carrier
        plot(time, carrier_velocity);
        % Make plot look pretty
        grid on;
        xlabel("Time [s]"); ylabel("Velocity [m/s]");
    end
end
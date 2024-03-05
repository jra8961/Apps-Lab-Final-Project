% MECE 301 - Engineering Applications Lab
% Final Project - Pneumatic Tube
% Written by Jacob, Jason, Jonah, and Tim

% Clear command window, variables, and figures
clc; clear; close all;

%% PHYSICAL CONSTANTS
% USE FORMAT : [estimate_value uncertainty]

% Room Temp and Atmospheric Pressure Air
room_temp = [22 1]; % C
p_atmos = [101.3 1] * 10^3; % Pa

% Dynamic viscosity of air
% SOURCE: https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html
air_viscosity = [1.822 0.0001822] *10^-5; %N*s/m2

% Tube Dimensions
tube_diameter = [1.75 0.012] * 0.0254; % in -> m
tube_length = [48 0.5] * 0.0254; % in -> m

% Carrier Dimensions
carrier_diameter = [1.735 0.002] * 0.0254; % in -> m
carrier_mass = [0.250 0.001]; % kg
min_carrier_length = 0.05; % m
max_carrier_length = 0.15; % m

%% RUN SIMULATION
% Number of points on the length vs. velocity graph
n_pts = 20;

% Setup arrays for storing carrier length and terminal velocity
terminal_velocity = zeros(n_pts,3);
carrier_length = linspace(min_carrier_length,max_carrier_length,n_pts); % kPa

% Store all physical parameters in an array to make it easy to pass into function
system_parameters = [room_temp(1), p_atmos(1), air_viscosity(1), tube_diameter(1), tube_length(1), ...
    carrier_diameter(1), carrier_mass(1)];

% Store all variations in system parameters for easy access during loop
system_variations = [room_temp(2), p_atmos(2), air_viscosity(2), tube_diameter(2), tube_length(2), ...
    carrier_diameter(2), carrier_mass(2)];

% Loop for each carrier length
for m = 1:n_pts
    % Calculate estimate exit velocity for that vacuum pressure
    terminal_velocity(m,2) = modelTube(carrier_length(m),system_parameters,false);
    
    % Initialize variables for storing variations caused by each parameter
    max_error = 0;
    min_error = 0;
    
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

%% PLOTTING
% Sample plot of carrier position, velocity, and acceleration
modelTube(max_carrier_length,system_parameters,true);

% Flip velocity to make the plot easier to read
terminal_velocity = -terminal_velocity;

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
ylim([0 0.2]);
legend("Location","northeast");

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
function [term_vel] = modelTube(carrier_len, sys_params, printing)

% GIVEN CONSTANTS - taken from input array
   % Gravitational acceleration
    g = 9.81; % m/s^2 

    room_T = sys_params(1);
    p_atm = sys_params(2);
    viscosity = sys_params(3);
    tube_diam = sys_params(4);
    tube_len = sys_params(5);
    carrier_diam = sys_params(6);
    carrier_mass = sys_params(7);

% CALCULATED CONSTANTS  
    % Gap between the carrier and tube
    air_gap = (tube_diam - carrier_diam) / 2;

    % Density of air at room temperature and atmospheric pressure
    room_density = calcAirDensity(p_atm, room_T);

% TIME DEPENDENT VARIABLES
    n_steps = 90000; % Number of steps in the loop
    delta_t = 0.0001;  % Amount of time between each step

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
    p_inside(1) = p_atm;
    inside_volume(1) = tube_len * pi()/4 * tube_diam^2;
    inside_mass(1) = inside_volume(1) * room_density;
    inside_density(1) = room_density;

% MAIN EULER CALCULATION LOOP
    n = 1;
    while n < n_steps
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
        accel(n) = -g + (f_pressure + f_shear) / carrier_mass;

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
        
        if(exited_mass > mass_displaced)
            ans = 1;
        end

        % How much air is left underneath the carrier?
        inside_mass(n+1) = inside_mass(n) - exited_mass;
        inside_volume(n+1) = (tube_len + position(n)) * pi()/4*tube_diam^2;
        inside_density(n+1) = inside_mass(n+1) / inside_volume(n+1);

        % Use ideal gas law to determine what the next pressure
        p_inside(n+1) = inside_volume(1)/inside_volume(n) * inside_mass(n)/inside_mass(1) * p_inside(1);

        % Increment step counter
        n = n + 1;
    end
    
% PLACEHOLDER !!!! CHANGE LATER
    term_vel = carrier_velocity(end);

% SAMPLE PLOT PRINTING
    if printing
        figure();
        
        % Plot position of carrier
        subplot(3,1,1);
        plot(time,position);
        % Make plot look pretty
        grid on;
        xlabel("Time [s]"); ylabel("Position [m]");
        title("POSITION");

        % Plot velocity of the carrier
        subplot(3,1,2);
        plot(time, carrier_velocity);
        % Make plot look pretty
        grid on;
        xlabel("Time [s]"); ylabel("Velocity [m/s]");
        title("VELOCITY");

        % Plot another thing - You decide!
        subplot(3,1,3);
        plot(time,accel);
        % Make plot look pretty
        grid on;
        xlabel("Time [s]"); ylabel("m/s2"); 
        title("ACCEL");
    end
end
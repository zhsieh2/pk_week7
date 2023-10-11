% Parameters
Vd = 3;  % Volume of distribution (in ml)
C0 = 100;  % Initial concentration (in mg/ml)
k0 = 15;  % Zero-order elimination rate (in mg/hour)
% Time parameters
time_injection = 0;  % Start of the initial IV injection
time_infusion = 24;  % Start of the IV infusion on the next day
infusion_duration = 6;  % Duration of the infusion (in hours)

% Define the time span for the simulation
tspan = [0, 30];  % 30 hours in total
% Define the function for the differential equation (zero-order kinetics)
ode = @(t, y) -k0;  % Constant elimination rate
% Solve the differential equation using ODE45
[t, concentrations] = ode45(ode, tspan, C0);

% Adjust concentrations during the IV infusion period
infusion_indices = t >= time_infusion & t <= time_infusion + infusion_duration;
infusion_start = find(t == time_infusion);
concentrations(infusion_indices)=concentrations(infusion_indices)+(t(infusion_indices)-time_infusion)*infusion_rate();

% Set concentrations to 0 if they go below 0 during the initial injection
concentrations(concentrations < 0) = 0;

% Plot the drug concentration over time
plot(t, concentrations);
xlabel('Time (hours)');
ylabel('Drug Concentration (mg/ml)');
title('Drug Concentration over Time');
grid on;

figure();
% Plot the log(drug concentration) over time
plot(t, log(concentrations));
xlabel('Time (hours)');
ylabel('log(Drug Concentration)');
title('log(Drug Concentration) over Time');
grid on;

% Function for the infusion rate
function rate = infusion_rate()
    rate = 100;  % ml/hour
end


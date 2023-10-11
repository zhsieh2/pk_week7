% Define the differential equation for drug concentration
% dC/dt = (Rate of Infusion - Rate of Elimination) / Volume
% We'll use zeroth-order elimination kinetics

% Constants
initialVolume = 3; % ml
infusionRate = 100; % ml/hour
infusionConcentration = 500; % mg/liter
eliminationRate = 15; % mg/hour

% Time span for the simulation
tspan_initial = [0 24]; % 0 to 24 hours for initial IV injection
tspan_infusion = [24 48]; % 24 to 48 hours for IV infusion

% Initial condition for IV injection (drug concentration at t=0)
C0_initial = 100; % Initial concentration after IV injection

% Initial condition for IV infusion (drug concentration at t=24 hours)
C0_infusion = 0; % No drug left at the start of IV infusion

% Define the differential equation as a function
ode = @(t, C) (infusionRate * infusionConcentration - eliminationRate) / initialVolume;

% Solve the differential equation for IV injection using ode45
[t_initial, C_initial] = ode45(ode, tspan_initial, C0_initial);

% Solve the differential equation for IV infusion using ode45
[t_infusion, C_infusion] = ode45(ode, tspan_infusion, C0_infusion);

% Combine the results for both IV injection and IV infusion
t_combined = [t_initial; t_infusion];
C_combined = [C_initial; C_infusion];

% Plot the drug concentration over time
plot(t_combined, C_combined, 'b', 'LineWidth', 2);
xlabel('Time (hours)');
ylabel('Drug Concentration (mg/ml)');
title('Drug Concentration Over Time');
grid on;
hold off;

figure();
% Plot the log(drug concentration) over time
plot(t_combined, log(C_combined), 'r', 'LineWidth', 2);
xlabel('Time (hours)');
ylabel('log(Drug Concentration)');
title('log(Drug Concentration) Over Time');
grid on;

% Calculate AUC for the IV infusion
infusionVolume = infusionRate * infusionConcentration * 6 / 1000; % Convert ml to liters
AUC_infusion = infusionVolume;

% Display the AUC value for the IV infusion
disp(['AUC for the IV infusion: ' num2str(AUC_infusion) ' mg']);


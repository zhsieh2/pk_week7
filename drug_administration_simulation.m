function drug_administration_simulation()
    % Parameters for the two-compartment model
    V1 = 100;  % Volume of distribution in central compartment (mL)
    V2 = 200;  % Volume of distribution in peripheral compartment (mL)
    CL = 10;   % Clearance (mL/min)
    
    % Absorption rate constants for PO, IM, and SC
    ka_PO = 0.2;  % Absorption rate constant for PO (1/min)
    ka_IM = 0.3;  % Absorption rate constant for IM (1/min)
    ka_SC = 0.1;  % Absorption rate constant for SC (1/min)
    
    % Distribution rate constants
    k12 = 0.05;  % Rate of drug flow from central to peripheral (1/min)
    k21 = 0.1;   % Rate of drug flow from peripheral to central (1/min)

    % Simulation parameters
    tspan = [0 60];  % Simulation time (minutes)

    % Initial drug concentrations in central and peripheral compartments for IV, PO, IM, and SC
    initial_IV = [100; 0];  % IV administration
    initial_PO = [0; 0];  % PO administration
    initial_IM = [0; 0];  % IM administration
    initial_SC = [0; 0];  % SC administration
    
    % Absorption durations proportional to the absorption constants
    absorption_duration_PO = 1 / ka_PO;  % Proportional to 1/ka_PO
    absorption_duration_IM = 1 / ka_IM;  % Proportional to 1/ka_IM
    absorption_duration_SC = 1 / ka_SC;  % Proportional to 1/ka_SC

    % Solve the differential equations using ode45 for each administration route
    [t_IV, drug_conc_IV] = ode45(@(t, y) ode_equations(t, y, V1, V2, CL, 0, 0, k12, k21, 0), tspan, initial_IV);
    [t_PO, drug_conc_PO] = ode45(@(t, y) ode_equations(t, y, V1, V2, CL, ka_PO, 0, k12, k21, absorption_duration_PO), tspan, initial_PO);
    [t_IM, drug_conc_IM] = ode45(@(t, y) ode_equations(t, y, V1, V2, CL, ka_IM, 0, k12, k21, absorption_duration_IM), tspan, initial_IM);
    [t_SC, drug_conc_SC] = ode45(@(t, y) ode_equations(t, y, V1, V2, CL, ka_SC, 0, k12, k21, absorption_duration_SC), tspan, initial_SC);

    % Plot drug concentrations for each administration route
    plot_drug_concentrations(t_IV, drug_conc_IV, 'IV Administration');
    plot_drug_concentrations(t_PO, drug_conc_PO, 'PO Administration');
    plot_drug_concentrations(t_IM, drug_conc_IM, 'IM Administration');
    plot_drug_concentrations(t_SC, drug_conc_SC, 'SC Administration');
end

function plot_drug_concentrations(t, drug_conc, title_str)
    % Plot drug concentrations and log(drug concentrations) over time
    figure;
    subplot(2, 1, 1);
    plot(t, drug_conc(:, 1), 'r-', t, drug_conc(:, 2), 'b-');
    xlabel('Time (min)');
    ylabel('Drug Concentration (ng/mL)');
    legend('Central Compartment', 'Peripheral Compartment');
    title(['Drug Concentration vs. Time (' title_str ')']);
    
    % Calculate and display AUC for drug concentrations
    auc = trapz(t, drug_conc(:, 1));
    disp(['AUC for ' title_str ' (Drug Concentrations): ' num2str(auc)]);
    
    subplot(2, 1, 2);
    plot(t, log(drug_conc(:, 1)), 'r-', t, log(drug_conc(:, 2)), 'b-');
    xlabel('Time (min)');
    ylabel('Log Drug Concentration (ng/mL)');
    legend('Log Central Compartment', 'Log Peripheral Compartment');
    title(['Log Drug Concentration vs. Time (' title_str ')']);
    
   % Calculate and display AUC for drug concentrations in the central compartment
    auc_central = trapz(t, drug_conc(:, 1));
    disp(['AUC for ' title_str ' (Central Compartment): ' num2str(auc_central)]);
    
    % Calculate and display AUC for drug concentrations in the peripheral compartment
    auc_peripheral = trapz(t, drug_conc(:, 2));
    disp(['AUC for ' title_str ' (Peripheral Compartment): ' num2str(auc_peripheral)]);
    
end

function dydt = ode_equations(t, y, V1, V2, CL, ka, kb, k12, k21, absorption_duration)
    % Differential equations for the two-compartment model with different administration routes
    C1 = y(1);  % Concentration in the central compartment (ng/mL)
    C2 = y(2);  % Concentration in the peripheral compartment (ng/mL)

    % Absorption (ka) for PO, IM, and SC within the absorption duration
    if absorption_duration > 0 && t <= absorption_duration
        if ka > 0
            dydt(1, 1) = ka * 100 - CL/V1 * C1 - k12 * C1 + k21 * C2;  % Rate of change in the central compartment
                                                                       % assume: drug concentration when uptaken = 100
            dydt(2, 1) = k12 * C1 - k21 * C2;  % Rate of change in the peripheral compartment
        else
            % Rate of change of drug concentrations in each compartment based on administration route
            dydt(1, 1) = -CL/V1 * C1 - k12 * C1 + k21 * C2;  % Rate of change in the central compartment
            dydt(2, 1) = k12 * C1 - k21 * C2;  % Rate of change in the peripheral compartment
        end
    else
        % No absorption beyond the absorption duration
        dydt(1, 1) = -CL/V1 * C1 - k12 * C1 + k21 * C2;  % Rate of change in the central compartment
        dydt(2, 1) = k12 * C1 - k21 * C2;  % Rate of change in the peripheral compartment
    end
end

%% ========================================================================
%  MAIN KINETIC MODELING WORKFLOW (CO2RR)
%  Global Optimization using MultiStart
% =========================================================================
clear; clc; close all;
tic;

%% ========================================================================
%  SECTION 1: Data Loading and Preparation
% =========================================================================
fprintf('Loading experimental data...\n');

% --- [USER INPUT REQUIRED] ---
% Specify the filenames for your experimental data. 
% Ensure these files are located in the same directory as this script.
file_currents = 'data_current_LSV.xlsx';       % File containing V, j_CO, j_H2
file_concentrations = 'data_concentrations.xlsx'; % File containing CO2, OH, pH
try
    % --- 1. Load Current Density Data ---
    % Assuming a simple table format. Modify the columns as needed based on your Excel structure.
    temp_table = readtable(file_currents);
    V_Matrix = temp_table{:, 1};           % Uncompensated Applied Potential (V)
    j_CO_Matrix = temp_table{:, 2} * 10;   % CO partial current density (unit conversion if needed)
    j_H2_Matrix = temp_table{:, 3} * 10;   % H2 partial current density (unit conversion if needed)

    % --- 2. Load Local Concentration Data ---
    C.CO2 = table2array(readtable(file_concentrations, 'Sheet', 'CO2'));
    C.OH  = table2array(readtable(file_concentrations, 'Sheet', 'OH'));
    C.pH  = table2array(readtable(file_concentrations, 'Sheet', 'pH'));

    % --- 3. Physical Constants and Corrections ---
    F = 96485;
    Rs = 9.2; % Average Series Resistance (Ohms) - Update this value based on your cell

    % --- 4. Process Data for Optimization ---
    % Calculate total current for IR compensation
    jtotal = j_CO_Matrix + j_H2_Matrix;
    
    % Correct Potential (IR compensation, Reference Electrode shift, and pH correction)
    % V = -(V_applied + iR_drop + E_reference + Nernst_pH_shift)
    V = -(V_Matrix + (jtotal) ./ 10000 .* Rs * 0.70 + 0.1971 + 0.059 .* C.pH); 

    % Convert current densities to experimental reaction rates (r)
    r.CO = j_CO_Matrix ./ (2 * F);
    r.H2 = j_H2_Matrix ./ (2 * F);

    fprintf('Data successfully loaded and formatted for optimization.\n');

catch ME
    fprintf('\nERROR: Failed to load data.\n');
    fprintf('Please ensure "%s" and "%s" exist in the current folder.\n', file_currents, file_concentrations);
    error(ME.message);
end

%% ========================================================================
%  SECTION 2: Global Optimization Setup
% =========================================================================
% Set num_runs = 2000 for production, reduced to 100 here for quick testing
num_runs = 2000; 

% Initialize parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool('local'); 
end

fprintf('Setting up global optimization problem for 12 parameters (8 k''s, 4 a''s)...\n');

% --- Define Problem Constraints and Bounds ---
k_LB = ones(1, 8) * -13;    % Lower bound for log10(k)
k_UB = ones(1, 8) * 8;      % Upper bound for log10(k)
a_LB = ones(1, 4) * 0.1;    % Lower bound for alpha
a_UB = ones(1, 4) * 0.9;    % Upper bound for alpha

LB = [k_LB, a_LB]';
UB = [k_UB, a_UB]';

% --- Initial Guess ---
initial_guess_k = [1.19e-06, 0.229, 1.55e-05, 0.00569, 0.00114, 4.169e-08, 0.417, 8.27e-08];
initial_guess_a = [0.5, 0.5, 0.5, 0.5];
X0 = [log10(initial_guess_k), initial_guess_a]';

% --- Objective Function Handle ---
objective_fun = @(params) ObjectiveFun_Log(V, C, r, params);

% --- Solver Options ---
options = optimoptions('lsqnonlin', ...
    'Algorithm','trust-region-reflective', ...
    'MaxIterations', 1e4, ...
    'MaxFunctionEvaluations', 5e5, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'Display', 'none');
% --- Create the Optimization Problem Structure ---
problem = createOptimProblem('lsqnonlin', ...
    'objective', objective_fun, ...
    'x0', X0, ...
    'lb', LB, ...
    'ub', UB, ...
    'options', options);

% --- MultiStart Configuration ---
ms = MultiStart('UseParallel', true, 'Display', 'iter');

% --- Generate Custom Start Points ---
fprintf('Generating a custom, centrally-focused set of %d start points...\n', num_runs);
num_params = length(X0);
custom_points = zeros(num_params, num_runs);
search_range = UB - LB;
stdev = search_range * 0.15; % 15% standard deviation

for i = 1:num_runs
    while true
        new_point = X0 + stdev .* randn(num_params, 1);
        if all(new_point >= LB) && all(new_point <= UB)
            custom_points(:, i) = new_point;
            break; 
        end
    end
end
start_points = CustomStartPointSet(custom_points');

%% ========================================================================
%  SECTION 3: Execution and Analysis
% =========================================================================
fprintf('Starting MultiStart Global Search...\n');
[X_solution_best, fval_best, exitflag_best, output, all_solutions] = run(ms, problem, start_points);

fprintf('\n=== Optimization Finished ===\n');
fprintf('Best Sum of Squared Residuals: %e\n', fval_best);

% Extract optimized parameters
log10_k_opt = X_solution_best(1:8);
a_opt = X_solution_best(9:12);
X_opt_k = 10.^log10_k_opt;

%% ========================================================================
%  SECTION 4: Simulate Coverages and Current Densities 
% =========================================================================
f = F / (8.314 * 298.15);
a1 = a_opt(1); a2 = a_opt(2); a4 = a_opt(3); a5 = a_opt(4);
aH2O = 1; aOH = C.OH; aCO2 = C.CO2;

k1f = X_opt_k(1).*exp(a1.*f.*V);
k1b = X_opt_k(2).*exp(-(1-a1).*f.*V);
k2f = X_opt_k(3).*exp(a2.*f.*V); k2b = X_opt_k(4).*exp(-(1-a2).*f.*V);
k3f = X_opt_k(5);  
k4f = X_opt_k(6).*exp(a4.*f.*V);
k4b = X_opt_k(7).*exp(-(1-a4).*f.*V);
k5f = X_opt_k(8).*exp(a5.*f.*V);

epsilon = 1e-25;
denom_H = k4b.*aOH + k5f.*aH2O + epsilon;
denom_CO_main = k2f.*k3f + k1b.*k3f.*aOH + k1b.*k2b.*(aOH).^2 + epsilon;

term_H = (k4f.*aH2O) ./ denom_H;
term_CO = (k1f.*aCO2.*aH2O.*(k2f + k3f + k2b.*aOH)) ./ denom_CO_main;

% Simulated Coverages
Results.cBare = 1./(term_H + term_CO + 1);
Results.cCOOH = ((k1f.*aCO2.*aH2O.*(k3f + k2b.*aOH)) ./ denom_CO_main) .* Results.cBare;
Results.cCO = ((k1f.*k2f.*aCO2.*aH2O) ./ denom_CO_main) .* Results.cBare;
Results.cH = (term_H) .* Results.cBare;

% Simulated Current Densities
Pred.jCO = 2 .* F .* (k3f .* Results.cCO);
Pred.jH2 = 2 .* F .* (k5f .* Results.cH .* aH2O);

% Save Results
filename = 'kinetics_optimization_results.mat'; 
save(filename, 'Results', 'Pred', 'X_opt_k', 'a_opt', 'fval_best');
fprintf('Results saved to %s\n', filename);

toc;
function Residuals = ObjectiveFun_Log(V, C, r, params)
% OBJECTIVEFUN_LOG Calculates the log-space residuals for CO2RR kinetic modeling.
% 
% Inputs:
%   V      - Applied potential (V vs. reference)
%   C      - Structure containing local concentrations (C.CO2, C.OH)
%   r      - Structure containing experimental reaction rates (r.CO, r.H2)
%   params - 12x1 vector containing 8 kinetic constants (log10_k) and 4 transfer coefficients (a)
%
% Output:
%   Residuals - Column vector of weighted residuals for CO and H2 fitting

    % --- 1. Define Weights ---
    w_CO = 5;
    w_H2 = 1; 
    
    % --- 2. Unpack Parameters (12 total) ---
    % First 8 parameters are log10(k), last 4 are alpha (transfer coefficients)
    log10_k = params(1:8);
    a_coeffs = params(9:12);

    vars = 10.^log10_k;
    a1 = a_coeffs(1);
    a2 = a_coeffs(2);
    a4 = a_coeffs(3);
    a5 = a_coeffs(4);
    
    % --- 3. Physical Constants ---
    F = 96485;                  % Faraday constant (C/mol)
    f = F / (8.314 * 298.15);   % f = F/RT
    aH2O = 1;                   % Activity of water
    aOH = C.OH;
    aCO2 = C.CO2;
    
    % Floor to prevent division-by-zero or log(0)
    epsilon = 1e-25;

    % --- 4. Calculate Rate Constants (k values) ---
    k1f = vars(1) .* exp(a1 .* f .* V); 
    k1b = vars(2) .* exp(-(1-a1) .* f .* V); 
    k2f = vars(3) .* exp(a2 .* f .* V); 
    k2b = vars(4) .* exp(-(1-a2) .* f .* V);
    k3f = vars(5);  
    k4f = vars(6) .* exp(a4 .* f .* V);
    k4b = vars(7) .* exp(-(1-a4) .* f .* V);
    k5f = vars(8) .* exp(a5 .* f .* V);
    
    % --- 5. Calculate Surface Coverages for Stability ---
    denom_H = k4b .* aOH + k5f .* aH2O + epsilon;
    denom_CO_main = k2f .* k3f + k1b .* k3f .* aOH + k1b .* k2b .* (aOH).^2 + epsilon;
    
    term_H = (k4f .* aH2O) ./ denom_H;
    term_CO = (k1f .* aCO2 .* aH2O .* (k2f + k3f + k2b .* aOH)) ./ denom_CO_main;

    c_bare = 1 ./ (term_H + term_CO + 1);
    c_CO = ((k1f .* k2f .* aCO2 .* aH2O) ./ denom_CO_main) .* c_bare;
    c_H = (term_H) .* c_bare;

    % --- 6. Calculate Final Model Rates ---
    r_CO_model = k3f .* c_CO;
    r_H2_model = k5f .* c_H .* aH2O;

    % --- 7. Calculate Residuals in Log-Space ---
    log_r_CO_model = log(max(r_CO_model, epsilon));
    log_r_CO_exp   = log(max(r.CO, epsilon));

    log_r_H2_model = log(max(r_H2_model, epsilon));
    log_r_H2_exp   = log(max(r.H2, epsilon));

    res_CO = w_CO * (log_r_CO_model - log_r_CO_exp);
    res_H2 = w_H2 * (log_r_H2_model - log_r_H2_exp);

    Residuals = [res_CO(:); res_H2(:)];

    % --- 8. Robustness Check ---
    if any(isnan(Residuals)) || any(isinf(Residuals))
        Residuals = ones(size(Residuals)) * 1e10; % Return a huge error penalty
    end
end
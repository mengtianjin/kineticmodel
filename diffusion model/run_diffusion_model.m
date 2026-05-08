function results = run_diffusion_model(config)
%RUN_DIFFUSION_MODEL  Solve the CO2RR diffusion-reaction model.
%
%   results = run_diffusion_model(config)
%
%   Required inputs (not included in the repository):
%     config.bulkFile    - Excel file containing bulk concentrations
%     config.currentFile - Excel file containing current density / voltage data
%
%   The script preserves the main algorithm needed to reproduce the local
%   concentration results.

    if nargin < 1 || isempty(config)
        config = default_config();
    end

    assert(isfile(config.bulkFile),   'Bulk concentration file not found: %s', config.bulkFile);
    assert(isfile(config.currentFile),'Current-density file not found: %s', config.currentFile);

    if ~isfolder(config.outputDir)
        mkdir(config.outputDir);
    end

    bulkData = readmatrix(config.bulkFile); % calculate the bulk concentration firtst
    p_vector = bulkData(:,1).';
    bulk_a_master = bulkData(:,3).' * 1e3;   % CO2
    bulk_c_master = bulkData(:,5).' * 1e3;   % HCO3-
    bulk_d_master = bulkData(:,6).' * 1e3;   % CO3^2-
    bulk_y_master = bulkData(:,8).' * 1e3;   % OH-

    if exist('sheetnames', 'file')
        sheetNames = sheetnames(config.currentFile);
    else
        [~, sheetNames] = xlsfinfo(config.currentFile);
        sheetNames = string(sheetNames);
    end
    numSheets = numel(sheetNames);

    [V_matrix, j_CO_matrix, j_HAC_matrix, j_H2_matrix] = read_current_workbook(config.currentFile, sheetNames);

    nRows = size(V_matrix, 1);
    nCols = numSheets;

    % Spatial and Time Mesh
    xmesh = linspace(0, config.L, config.nx);
    tspan = linspace(0, config.tEnd, config.nt);

    % ODE solver options for pdepe to enforce non-negativity
    % This applies to all 4 components (a, d, c, y)
    pdeOptions = odeset('NonNegative', 1:4, 'RelTol', config.relTol, 'AbsTol', config.absTol);

    %% Main Simulation Loop
    allRunsData = cell(nRows, nCols);
    A_end = zeros(nRows, nCols);
    C_end = zeros(nRows, nCols);
    D_end = zeros(nRows, nCols);
    Y_end = zeros(nRows, nCols);

    for rowIdx = 1:nRows
        for colIdx = 1:nCols
            bulk_a = max(0, bulk_a_master(min(colIdx,end)));
            bulk_c = max(0, bulk_c_master(min(colIdx,end)));
            bulk_d = max(0, bulk_d_master(min(colIdx,end)));
            bulk_y = max(0, bulk_y_master(min(colIdx,end)));

            j_CO  = j_CO_matrix(rowIdx, colIdx);
            j_H2  = j_H2_matrix(rowIdx, colIdx);
            j_HAC = j_HAC_matrix(rowIdx, colIdx);

            flux_a = -(j_HAC / (2 * config.F)+j_CO / (2 * config.F));
            flux_y = -(j_CO / (2 * config.F)+j_HAC / (2 * config.F) + j_H2 / config.F);
         
            % Solve PDE System with NonNegative option
            sol = pdepe(0, ...
                @(x,t,u,DuDx) pde_system(x, t, u, DuDx, config), ...
                @(x) initial_conditions(x, bulk_a, bulk_d, bulk_c, bulk_y), ...
                @(xl,ul,xr,ur,t) boundary_conditions(xl, ul, xr, ur, t, bulk_a, bulk_d, bulk_c, bulk_y, flux_a, flux_y), ...
                xmesh, tspan, pdeOptions);

            runData.rowIdx = rowIdx;
            runData.colIdx = colIdx;
            runData.voltage = V_matrix(rowIdx, colIdx);
            runData.pressure = p_vector(min(colIdx, numel(p_vector)));
            runData.xmesh = xmesh;
            runData.tspan = tspan;
            runData.a = sol(:,:,1);
            runData.d = sol(:,:,2);
            runData.c = sol(:,:,3);
            runData.y = sol(:,:,4);

            allRunsData{rowIdx, colIdx} = runData;

            A_end(rowIdx, colIdx) = sol(end, end, 1);
            D_end(rowIdx, colIdx) = sol(end, end, 2);
            C_end(rowIdx, colIdx) = sol(end, end, 3);
            Y_end(rowIdx, colIdx) = sol(end, end, 4);
        end
    end

    results.allRunsData = allRunsData;
    results.V_matrix = V_matrix;
    results.p_vector = p_vector;
    results.A_end = A_end;
    results.C_end = C_end;
    results.D_end = D_end;
    results.Y_end = Y_end;
    results.pH_end = 14 + log10(max(Y_end, eps) * 1e-3);
    results.xmesh = xmesh;
    results.tspan = tspan;

    save(fullfile(config.outputDir, config.outputFile), 'results');

    fprintf('Model runs completed. Results saved to %s\n', fullfile(config.outputDir, config.outputFile));
end

function config = default_config()
    rootDir = pwd;

    config.bulkFile = fullfile(rootDir, 'data', 'bulk.xlsx');
    config.currentFile = fullfile(rootDir, 'data', 'partial_current.xlsx');
    config.outputDir = fullfile(rootDir, 'results');
    config.outputFile = 'diffusion_model_results.mat';

    config.F = 96485;
    config.L = 5e-5;     % diffusion layer thickness 
    config.nx = 100;
    config.tEnd = 20;
    config.nt = 50;
    config.relTol = 1e-5;
    config.absTol = 1e-8;

    % Kinetic constants
    config.k7f = 5.93;
    config.k7b = 1.34e-4;
    config.k8f = 1e5;
    config.k8b = 2.15e4;

    % Diffusion coefficients
    config.Da = 1.91e-9;
    config.Dd = 9.23e-10;
    config.Dc = 1.19e-9;
    config.Dy = 5.27e-9;
end

function [V_matrix, j_CO_matrix, j_HAC_matrix, j_H2_matrix] = read_current_workbook(fileName, sheetNames)
    numSheets = numel(sheetNames);
    firstSheet = readtable(fileName, 'Sheet', sheetNames{1});
    nRows = height(firstSheet);

    V_matrix = zeros(nRows, numSheets);
    j_CO_matrix = zeros(nRows, numSheets);
    j_HAC_matrix = zeros(nRows, numSheets);
    j_H2_matrix = zeros(nRows, numSheets);

    for i = 1:numSheets
        temp = readtable(fileName, 'Sheet', sheetNames{i});
        tempArray = table2array(temp);

        n = min(nRows, size(tempArray,1));
        V_matrix(1:n, i) = tempArray(1:n, 1);
        j_CO_matrix(1:n, i) = tempArray(1:n, 4) * 10;
        j_HAC_matrix(1:n, i) = tempArray(1:n, 2) * 10;
        j_H2_matrix(1:n, i) = tempArray(1:n, 3) * 10;
    end
end
%% Helper Functions (Unchanged)
function [c_coeff, f_flux, s_source] = pde_system(~, ~, u, DuDx, cfg)
    c_coeff = [1; 1; 1; 1];
    f_flux = [cfg.Da * DuDx(1); cfg.Dd * DuDx(2); cfg.Dc * DuDx(3); cfg.Dy * DuDx(4)];

    a = max(0, u(1));
    d = max(0, u(2));
    c = max(0, u(3));
    y = max(0, u(4));

    s_a = cfg.k7b * c - cfg.k7f * a * y;
    s_d = cfg.k8f * c * y - cfg.k8b * d;
    s_c = cfg.k7f * a * y + cfg.k8b * d - cfg.k7b * c - cfg.k8f * c * y;
    s_y = cfg.k7b * c + cfg.k8b * d - cfg.k7f * a * y - cfg.k8f * c * y;

    s_source = [s_a; s_d; s_c; s_y];
end

function u0 = initial_conditions(~, bulk_a, bulk_d, bulk_c, bulk_y)
    u0 = [max(0, bulk_a); max(0, bulk_d); max(0, bulk_c); max(0, bulk_y)];
end

function [pl, ql, pr, qr] = boundary_conditions(~, ul, ~, ~, ~, bulk_a, bulk_d, bulk_c, bulk_y, flux_a, flux_y)
    % Left boundary: fixed bulk concentrations
    pl = [ul(1) - max(0, bulk_a); ...
          ul(2) - max(0, bulk_d); ...
          ul(3) - max(0, bulk_c); ...
          ul(4) - max(0, bulk_y)];
    ql = [0; 0; 0; 0];

    % Right boundary: flux conditions for species linked to electrochemistry
    pr = [-flux_a; 0; 0; -flux_y];
    qr = [1; 1; 1; 1];
end

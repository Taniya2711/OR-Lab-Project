function analyze_case14()
    % Ensure MATPOWER is in the path
    if ~exist('case14.m', 'file')
        error('MATPOWER not found. Add it to path: addpath(genpath(''matpower_folder''))');
    end

    % Load case data
    mpc = case14;
    
    % =============================================
    % 1. Extract Bus Data (P, Q, V_ref)
    % =============================================
    bus_data = [
        mpc.bus(:, 1), ...    % Bus number
        mpc.bus(:, 3), ...    % P_load (MW)
        mpc.bus(:, 4), ...    % Q_load (MVAr)
        mpc.bus(:, 8)         % V_ref (p.u.)
    ];
    bus_table = array2table(bus_data, ...
        'VariableNames', {'Bus', 'P_Load_MW', 'Q_Load_MVAr', 'V_ref_pu'});
    
    % =============================================
    % 2. Extract Branch Data (G, B, S_max)
    % =============================================
    % Calculate G and B from R and X
    R = mpc.branch(:, 3);     % Resistance (p.u.)
    X = mpc.branch(:, 4);     % Reactance (p.u.)
    G = R ./ (R.^2 + X.^2);   % Conductance (p.u.)
    B = -X ./ (R.^2 + X.^2);  % Susceptance (p.u.)
    S_max = mpc.branch(:, 6); % Thermal limit (MVA)
    
    branch_data = [
        mpc.branch(:, 1:2), ... % From/To buses
        G, B, S_max
    ];
    branch_table = array2table(branch_data, ...
        'VariableNames', {'From_Bus', 'To_Bus', 'G_pu', 'B_pu', 'S_max_MVA'});
    
    % =============================================
    % 3. Run Power Flow Analysis
    % =============================================
    results = runpf(mpc);
    
    % Extract power flow results
    voltage_magnitudes = results.bus(:, 8);     % Voltage (p.u.)
    line_flows = results.branch(:, 14);         % Real power flow (MW)
    
    % =============================================
    % 4. Display Results
    % =============================================
    disp('=== IEEE 14-Bus System Data ===');
    disp('1. Bus Data (Loads & Voltage):');
    disp(bus_table);
    
    disp('2. Branch Data (Admittance & Limits):');
    disp(branch_table);
    
    disp('3. Power Flow Results:');
    fprintf('Voltages: %s\n', mat2str(voltage_magnitudes, 3));
    fprintf('Line Flows (MW): %s\n', mat2str(line_flows, 3));
    
    % =============================================
    % 5. Export to Excel
    % =============================================
    writetable(bus_table, 'bus_data.xlsx');
    writetable(branch_table, 'branch_data.xlsx');
    disp('Data exported to bus_data.xlsx and branch_data.xlsx');
end

% Run the analysis
analyze_case14();
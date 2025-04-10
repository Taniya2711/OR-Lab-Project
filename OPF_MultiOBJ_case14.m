function analyze_opf_case14()
    % Load MATPOWER's 14-bus case and define constants
    mpc = case14;
    define_constants;
    
    % Ensure generator cost data exists (5 generators for case14)
    if size(mpc.gencost, 2) < 7
        mpc.gencost = [
            2, 0, 0, 2, 0.02, 10, 100;
            2, 0, 0, 2, 0.04, 8, 80;
            2, 0, 0, 2, 0.03, 12, 120;
            2, 0, 0, 2, 0.05, 15, 150;
            2, 0, 0, 2, 0.025, 9, 90;
        ];
    end
    % Set branch limits (RATE_A)
    mpc.branch(:, RATE_A) = 100;  % 100 MVA limit
    
    %% Extract System Parameters
    baseMVA = mpc.baseMVA;
    nb = size(mpc.bus, 1);      % Number of buses (14 for case14)
    ng = size(mpc.gen, 1);      % Number of generators (5 for case14)
    nl = size(mpc.branch, 1);   % Number of branches
    slack_bus = find(mpc.bus(:, BUS_TYPE) == 3);
    
    % Admittance matrix
    [Ybus, ~, ~] = makeYbus(baseMVA, mpc.bus, mpc.branch);
    G = real(Ybus);
    B = imag(Ybus);
    % Generator data
    Pmin = mpc.gen(:, PMIN) / baseMVA;
    Pmax = mpc.gen(:, PMAX) / baseMVA;
    Qmin = mpc.gen(:, QMIN) / baseMVA;
    Qmax = mpc.gen(:, QMAX) / baseMVA;
    if size(mpc.gencost, 2) >= 7
        a = mpc.gencost(:, 5);
        b = mpc.gencost(:, 6);
        c = mpc.gencost(:, 7);
    else
        a = 0.01 * ones(ng, 1);
        b = 10 * ones(ng, 1);
        c = 100 * ones(ng, 1);
    end
    Vmin = mpc.bus(:, VMIN);
    Vmax = mpc.bus(:, VMAX);
    Pd = mpc.bus(:, PD) / baseMVA;
    Qd = mpc.bus(:, QD) / baseMVA;
    %% Optimization Setup
    % Weighting factors for multi-objective (generation cost and voltage profile)
    w1 = 0.7;  % Weight on generation cost
    w2 = 0.3;  % Weight on voltage deviations
    Vref = 1.0;
    % Initial guess for decision variables: [V; theta; Pg; Qg]
    V0 = ones(nb, 1);
    theta0 = zeros(nb, 1);
    Pg0 = zeros(ng, 1);
    Qg0 = zeros(ng, 1);
    % Map generator bus indices & set mid-range generation values
    gen_buses = mpc.gen(:, 1);
    for i = 1:ng
        Pg0(i) = (Pmax(i) + Pmin(i)) / 2;
        Qg0(i) = (Qmax(i) + Qmin(i)) / 2;
    end
    x0 = [V0; theta0; Pg0; Qg0];
    % Define variable bounds
    lb_V = Vmin;
    ub_V = Vmax;
    lb_theta = -pi * ones(nb, 1);
    ub_theta = pi * ones(nb, 1);
    lb_theta(slack_bus) = 0;
    ub_theta(slack_bus) = 0;
    lb_Pg = Pmin;
    ub_Pg = Pmax;
    lb_Qg = Qmin;
    ub_Qg = Qmax;
    lb = [lb_V; lb_theta; lb_Pg; lb_Qg];
    ub = [ub_V; ub_theta; ub_Pg; ub_Qg];
    %% Define Objective and Constraint Functions
    objfun = @(x) opf_objective(x, nb, ng, gen_buses, a, b, c, w1, w2, Vref);
    nonlcon = @(x) opf_constraints(x, nb, ng, nl, gen_buses, G, B, Pd, Qd, slack_bus, mpc.branch(:, RATE_A)/baseMVA);
    %% Optimization Options and Solve with fmincon
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
        'Display', 'iter', 'MaxFunctionEvaluations', 10000, ...
        'MaxIterations', 1000, 'OptimalityTolerance', 1e-6, 'ConstraintTolerance', 1e-6);
    [x_opt, fval, exitflag, output] = fmincon(objfun, x0, [], [], [], [], lb, ub, nonlcon, options);
    %% Extract Results
    V_opt = x_opt(1:nb);
    theta_opt = x_opt(nb+1:2*nb);
    Pg_opt = x_opt(2*nb+1:2*nb+ng);
    Qg_opt = x_opt(2*nb+ng+1:end);
    %% Display Numerical Results
    disp('Optimal Bus Voltages and Angles:');
    for i = 1:nb
        fprintf('Bus %2d: V = %.4f p.u., theta = %.4f rad (%.2f°)\n', ...
            i, V_opt(i), theta_opt(i), theta_opt(i)*180/pi);
    end
    disp('Optimal Generator Outputs:');
    for i = 1:ng
        fprintf('Generator at Bus %2d: Pg = %.4f p.u. (%.2f MW), Qg = %.4f p.u. (%.2f MVAr)\n', ...
            gen_buses(i), Pg_opt(i), Pg_opt(i)*baseMVA, Qg_opt(i), Qg_opt(i)*baseMVA);
    end
    total_cost = sum(a .* Pg_opt.^2 + b .* Pg_opt + c);
    fprintf('Total Generation Cost: %.4f\n', total_cost);
    voltage_metric = sum((V_opt - Vref).^2);
    fprintf('Voltage Profile Metric: %.6f\n', voltage_metric);
    % losses = calculateLosses(V_opt, theta_opt, mpc.branch, baseMVA);
    % fprintf('Total System Losses: %.4f MW\n', losses);
    %% Enhanced Visual Outputs
    figure('Name', 'Optimal OPF Results - IEEE 14-Bus', 'Color', [1 1 1]);
    t = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    % Voltage Profile Plot
    nexttile;
    hBarVolt = bar(V_opt, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1.2);
    % Apply gradient coloring based on voltage levels
    for k = 1:length(V_opt)
        if V_opt(k) < 0.95
            hBarVolt.CData(k,:) = [0.8 0.2 0.2];  % Red tones for low voltages
        elseif V_opt(k) > 1.05
            hBarVolt.CData(k,:) = [0.2 0.2 0.8];  % Blue tones for high voltages
        else
            hBarVolt.CData(k,:) = [0.2 0.6 0.2];  % Green for acceptable voltage
        end
    end
    grid on; box on;
    xlabel('Bus Number');
    ylabel('Voltage (p.u.)');
    title('Optimal Bus Voltage Profile');
    ylim([min(V_opt)-0.05, max(V_opt)+0.05]);
    xticks(1:length(V_opt));
    % Annotate each bar with its voltage value
    for k = 1:length(V_opt)
        text(k, V_opt(k)+0.005, sprintf('%.3f', V_opt(k)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
    end
    % Voltage Angle Plot
    nexttile;
    hBarAngle = bar(theta_opt*180/pi, 'FaceColor', [0.6 0.4 0.8], 'EdgeColor', 'k', 'LineWidth', 1.2);
    grid on; box on;
    xlabel('Bus Number');
    ylabel('Angle (°)');
    title('Optimal Bus Voltage Angles');
    xticks(1:length(theta_opt));
    % Annotate each angle value (in degrees)
    for k = 1:length(theta_opt)
        text(k, (theta_opt(k)*180/pi)+0.5, sprintf('%.1f°', theta_opt(k)*180/pi), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
    end
    % Add an overall title for the tiled layout
    title(t, 'Enhanced Optimal OPF Results for IEEE 14-Bus System', 'FontSize', 14, 'FontWeight', 'bold');
end
%% Objective Function
function f = opf_objective(x, nb, ng, gen_buses, a, b, c, w1, w2, Vref)
    % Extract voltage and generator power values
    V = x(1:nb);
    Pg = x(2*nb+1:2*nb+ng);
    % Generation cost (quadratic model)
    gen_cost = sum(a .* Pg.^2 + b .* Pg + c);
    % Voltage deviation cost
    voltage_term = sum((V - Vref).^2);
    % Weighted combined objective
    f = w1 * gen_cost + w2 * voltage_term;
end
%% Nonlinear Constraints
function [c, ceq] = opf_constraints(x, nb, ng, nl, gen_buses, G, B, Pd, Qd, slack_bus, Smax)
    % Extract decision variables
    V = x(1:nb);
    theta = x(nb+1:2*nb);
    Pg = x(2*nb+1:2*nb+ng);
    Qg = x(2*nb+ng+1:2*nb+2*ng);
    % Map generator outputs to their buses
    Pg_bus = zeros(nb, 1);
    Qg_bus = zeros(nb, 1);
    for i = 1:ng
        Pg_bus(gen_buses(i)) = Pg(i);
        Qg_bus(gen_buses(i)) = Qg(i);
    end
    % Power balance equations for each bus
    Pbalance = zeros(nb, 1);
    Qbalance = zeros(nb, 1);
    for i = 1:nb
        P_inj = 0;
        Q_inj = 0;
        for j = 1:nb
            angle_diff = theta(i) - theta(j);
            P_inj = P_inj + V(i) * V(j) * (G(i,j) * cos(angle_diff) + B(i,j) * sin(angle_diff));
            Q_inj = Q_inj + V(i) * V(j) * (G(i,j) * sin(angle_diff) - B(i,j) * cos(angle_diff));
        end
        Pbalance(i) = Pg_bus(i) - Pd(i) - P_inj;
        Qbalance(i) = Qg_bus(i) - Qd(i) - Q_inj;
    end
    % Line flow constraints
    line_flow = [];
    if ~isempty(Smax) && Smax(1) ~= 999
        line_flow = computeLineFlows(V, theta, nb, nl, G, B, Smax);
    end
    c = line_flow;
    ceq = [Pbalance; Qbalance];
end
%% Compute Line Flows
function c = computeLineFlows(V, theta, nb, nl, G, B, Smax)
    c = [];
    line_count = 0;
    for i = 1:nb
        for j = i+1:nb
            if abs(G(i,j)) > 1e-6 || abs(B(i,j)) > 1e-6
                line_count = line_count + 1;
                if line_count <= nl
                    angle_diff = theta(i) - theta(j);
                    Pij = V(i)^2 * G(i,j) - V(i)*V(j) * (G(i,j)*cos(angle_diff) + B(i,j)*sin(angle_diff));
                    Qij = -V(i)^2 * B(i,j) - V(i)*V(j) * (G(i,j)*sin(angle_diff) - B(i,j)*cos(angle_diff));
                    Sij = sqrt(Pij^2 + Qij^2);
                    c = [c; Sij - Smax(line_count)];
                end
            end
        end
    end
end

%% Calculate Total System Losses
function losses = calculateLosses(V, theta, branch_data, baseMVA)
    losses = 0;
    for k = 1:size(branch_data, 1)
        from_bus = branch_data(k, 1);
        to_bus = branch_data(k, 2);
        r = branch_data(k, 3);
        x = branch_data(k, 4);
        angle_diff = theta(from_bus) - theta(to_bus);
        I_squared = (V(from_bus)^2 + V(to_bus)^2 - 2*V(from_bus)*V(to_bus)*cos(angle_diff)) / (r^2 + x^2);
        losses = losses + r * I_squared;
    end
    losses = losses * baseMVA;
end

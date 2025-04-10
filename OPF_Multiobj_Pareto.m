function analyze_opf_case14_pareto()
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
    %% Generate Pareto Front using epsilon-constraint method
    % Define number of points on the Pareto front
    num_points = 15;
    % Find the range of each objective
    % First, optimize for cost only
    objfun_cost = @(x) opf_cost_objective(x, nb, ng, gen_buses, a, b, c);
    nonlcon = @(x) opf_constraints(x, nb, ng, nl, gen_buses, G, B, Pd, Qd, slack_bus, mpc.branch(:, RATE_A)/baseMVA);
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
        'Display', 'iter', 'MaxFunctionEvaluations', 10000, ...
        'MaxIterations', 1000, 'OptimalityTolerance', 1e-6, 'ConstraintTolerance', 1e-6);
    [x_cost, fval_cost, ~, ~] = fmincon(objfun_cost, x0, [], [], [], [], lb, ub, nonlcon, options);
    % Extract cost-optimal results
    V_cost = x_cost(1:nb);
    voltage_dev_at_cost_opt = sum((V_cost - Vref).^2);
    fprintf('Min Cost Solution: Cost = %.4f, Voltage Deviation = %.6f\n', fval_cost, voltage_dev_at_cost_opt);
    % Then, optimize for voltage profile only
    objfun_voltage = @(x) opf_voltage_objective(x, nb, Vref);
    [x_volt, fval_volt, ~, ~] = fmincon(objfun_voltage, x0, [], [], [], [], lb, ub, nonlcon, options);
    % Extract voltage-optimal results
    cost_at_volt_opt = opf_cost_objective(x_volt, nb, ng, gen_buses, a, b, c);
    fprintf('Min Voltage Deviation Solution: Cost = %.4f, Voltage Deviation = %.6f\n', cost_at_volt_opt, fval_volt);
    % Define the range of epsilon values
    if voltage_dev_at_cost_opt > fval_volt
        epsilon_values = linspace(fval_volt, voltage_dev_at_cost_opt, num_points);
    else
        epsilon_values = linspace(fval_volt, 2*voltage_dev_at_cost_opt, num_points);
    end
    % Create arrays to store Pareto front results
    pareto_cost = zeros(num_points, 1);
    pareto_voltage_dev = zeros(num_points, 1);
    pareto_solutions = cell(num_points, 1);
    % Store cost-optimal solution
    pareto_cost(1) = fval_cost;
    pareto_voltage_dev(1) = voltage_dev_at_cost_opt;
    pareto_solutions{1} = x_cost;
    % Store voltage-optimal solution
    pareto_cost(num_points) = cost_at_volt_opt;
    pareto_voltage_dev(num_points) = fval_volt;
    pareto_solutions{num_points} = x_volt;
    % Solve for Pareto-optimal solutions using epsilon constraint method
    for i = 2:num_points-1
        epsilon = epsilon_values(i);
        fprintf('\nSolving for epsilon = %.6f (%d/%d)\n', epsilon, i, num_points);
        % Create constraint function with voltage deviation limit
        nonlcon_epsilon = @(x) opf_constraints_with_epsilon(x, nb, ng, nl, gen_buses, G, B, Pd, Qd, slack_bus, mpc.branch(:, RATE_A)/baseMVA, Vref, epsilon);
        % Solve with cost as objective and voltage deviation as constraint
        [x_pareto, fval_pareto, exitflag, ~] = fmincon(objfun_cost, x0, [], [], [], [], lb, ub, nonlcon_epsilon, options);
        if exitflag > 0
            % Store successful solution
            pareto_cost(i) = fval_pareto;
            pareto_voltage_dev(i) = opf_voltage_objective(x_pareto, nb, Vref);
            pareto_solutions{i} = x_pareto;
            fprintf('Found solution: Cost = %.4f, Voltage Deviation = %.6f\n', pareto_cost(i), pareto_voltage_dev(i));
        else
            fprintf('No solution found for epsilon = %.6f\n', epsilon);
            % Use previous solution
            pareto_cost(i) = pareto_cost(i-1);
            pareto_voltage_dev(i) = pareto_voltage_dev(i-1);
            pareto_solutions{i} = pareto_solutions{i-1};
        end
    end
    % Filter invalid or identical solutions
    idx = 1;
    valid_cost = pareto_cost(idx);
    valid_volt = pareto_voltage_dev(idx);
    valid_solutions = {pareto_solutions{idx}};
    for i = 2:num_points
        if ~isnan(pareto_cost(i)) && ~isnan(pareto_voltage_dev(i)) && ...
           (abs(pareto_cost(i) - valid_cost(end)) > 1e-4 || abs(pareto_voltage_dev(i) - valid_volt(end)) > 1e-6)
            valid_cost = [valid_cost; pareto_cost(i)];
            valid_volt = [valid_volt; pareto_voltage_dev(i)];
            valid_solutions = [valid_solutions; pareto_solutions(i)];
        end
    end
    %% Display and Visualize Pareto Front
    % Plot Pareto front
    figure('Name', 'Pareto Front - OPF Multi-Objective', 'Color', [1 1 1]);
    scatter(valid_cost, valid_volt, 100, 'filled', 'MarkerFaceColor', [0.3 0.6 0.9], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
    plot(valid_cost, valid_volt, 'k-', 'LineWidth', 1.2);
    xlabel('Generation Cost');
    ylabel('Voltage Deviation');
    title('Pareto Front: Generation Cost vs. Voltage Deviation');
    grid on; box on;
    text(valid_cost(1), valid_volt(1), '  Cost-optimal', 'FontWeight', 'bold');
    text(valid_cost(end), valid_volt(end), '  Voltage-optimal', 'FontWeight', 'bold');
    % Extract the solution for detailed analysis (mid-point of Pareto front as an example)
    mid_idx = ceil(length(valid_cost) / 2);
    x_selected = valid_solutions{mid_idx};
    % Extract values from selected solution
    V_opt = x_selected(1:nb);
    theta_opt = x_selected(nb+1:2*nb);
    Pg_opt = x_selected(2*nb+1:2*nb+ng);
    Qg_opt = x_selected(2*nb+ng+1:end);
    % Extract weights from the Pareto front
    extractWeightsFromPareto(valid_cost, valid_volt, valid_solutions);
    %% Display Numerical Results for selected solution
    fprintf('\n\nDetailed Results for Pareto Solution:\n');
    fprintf('Generation Cost: %.4f, Voltage Deviation: %.6f\n', valid_cost(mid_idx), valid_volt(mid_idx));
    fprintf('\nOptimal Bus Voltages and Angles:\n');
    for i = 1:nb
        fprintf('Bus %2d: V = %.4f p.u., theta = %.4f rad (%.2f°)\n', ...
            i, V_opt(i), theta_opt(i), theta_opt(i)*180/pi);
    end
    fprintf('\nOptimal Generator Outputs:\n');
    for i = 1:ng
        fprintf('Generator at Bus %2d: Pg = %.4f p.u. (%.2f MW), Qg = %.4f p.u. (%.2f MVAr)\n', ...
            gen_buses(i), Pg_opt(i), Pg_opt(i)*baseMVA, Qg_opt(i), Qg_opt(i)*baseMVA);
    end
    %% Enhanced Visual Outputs for Selected Solution
    figure('Name', 'Selected Pareto Solution - IEEE 14-Bus', 'Color', [1 1 1]);
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
    title(t, 'Selected Pareto-Optimal Solution for IEEE 14-Bus System', 'FontSize', 14, 'FontWeight', 'bold');
    %% Visualize trade-offs across Pareto front
    % Compare voltage profiles and generator outputs across Pareto front
    figure('Name', 'Pareto Solutions Comparison', 'Color', [1 1 1]);
    t2 = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    % Select solutions to compare (first, middle, last)
    compare_idx = [1, mid_idx, length(valid_cost)];
    solution_names = {'Cost-optimal', 'Mid-point', 'Voltage-optimal'};
    % Voltage profiles comparison
    nexttile;
    hold on;
    colors = [0.8 0.2 0.2; 0.2 0.6 0.2; 0.2 0.2 0.8]; % Red, Green, Blue
    markers = {'o-', 's-', 'diamond-'};
    for i = 1:length(compare_idx)
        idx = compare_idx(i);
        x_sol = valid_solutions{idx};
        V_sol = x_sol(1:nb);
        plot(1:nb, V_sol, markers{i}, 'LineWidth', 2, 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    end
    grid on; box on;
    xlabel('Bus Number');
    ylabel('Voltage (p.u.)');
    title('Voltage Profiles Across Pareto Front');
    xticks(1:nb);
    legend(solution_names, 'Location', 'best');
    ylim([min(Vmin)-0.02, max(Vmax)+0.02]);
    % Generator outputs comparison
    nexttile;
    gen_data = zeros(ng, length(compare_idx));
    for i = 1:length(compare_idx)
        idx = compare_idx(i);
        x_sol = valid_solutions{idx};
        Pg_sol = x_sol(2*nb+1:2*nb+ng);
        gen_data(:,i) = Pg_sol;
    end
    bar_h = bar(gen_data);
    for i = 1:length(compare_idx)
        bar_h(i).FaceColor = colors(i,:);
    end
    grid on; box on;
    xlabel('Generator Number');
    ylabel('Active Power (p.u.)');
    title('Generator Outputs Across Pareto Front');
    xticks(1:ng);
    xticklabels(arrayfun(@(x) sprintf('Gen %d', x), 1:ng, 'UniformOutput', false));
    legend(solution_names, 'Location', 'best');
    title(t2, 'Comparison of Different Solutions on the Pareto Front', 'FontSize', 14, 'FontWeight', 'bold');
end
%% Separate Objective Functions
function f = opf_cost_objective(x, nb, ng, gen_buses, a, b, c)
    % Cost-only objective
    Pg = x(2*nb+1:2*nb+ng);
    f = sum(a .* Pg.^2 + b .* Pg + c);
end
function f = opf_voltage_objective(x, nb, Vref)
    % Voltage deviation objective
    V = x(1:nb);
    f = sum((V - Vref).^2);
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
%% Constraints with epsilon constraint on voltage deviation
function [c, ceq] = opf_constraints_with_epsilon(x, nb, ng, nl, gen_buses, G, B, Pd, Qd, slack_bus, Smax, Vref, epsilon)
    % Regular constraints
    [c_reg, ceq] = opf_constraints(x, nb, ng, nl, gen_buses, G, B, Pd, Qd, slack_bus, Smax);
    % Epsilon constraint on voltage deviation
    V = x(1:nb);
    voltage_dev = sum((V - Vref).^2);
    c_eps = voltage_dev - epsilon;
    % Combine constraints
    c = [c_reg; c_eps];
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
%% Extract Weights from Pareto Front
function extractWeightsFromPareto(valid_cost, valid_volt, valid_solutions)
    fprintf('\n----- Weights Analysis from Pareto Front -----\n');
    fprintf('Point | Cost Value | Voltage Dev | w1 (Cost) | w2 (Voltage) | w1/w2 Ratio\n');
    fprintf('------------------------------------------------------------------\n');
    % Normalize the objectives to similar scales for fair comparison
    max_cost = max(valid_cost);
    min_cost = min(valid_cost);
    max_volt = max(valid_volt);
    min_volt = min(valid_volt);
    norm_cost = (valid_cost - min_cost) / (max_cost - min_cost);
    norm_volt = (valid_volt - min_volt) / (max_volt - min_volt);
    % Calculate weights for each point based on the local gradient
    for i = 1:length(valid_cost)
        if i == 1
            % First point - use forward difference
            dc = norm_cost(i+1) - norm_cost(i);
            dv = norm_volt(i+1) - norm_volt(i);
        elseif i == length(valid_cost)
            % Last point - use backward difference
            dc = norm_cost(i) - norm_cost(i-1);
            dv = norm_volt(i) - norm_volt(i-1);
        else
            % Interior points - use central difference
            dc = norm_cost(i+1) - norm_cost(i-1);
            dv = norm_volt(i+1) - norm_volt(i-1);
        end
        % Calculate gradient (slope of tangent line)
        if abs(dc) < 1e-10
            w1 = 0;
            w2 = 1;
        elseif abs(dv) < 1e-10
            w1 = 1;
            w2 = 0;
        else
            % Gradient of the Pareto front at this point
            gradient = -dv/dc;
            % Convert gradient to weights (w1 for cost, w2 for voltage
            w1 = 1 / (1 + gradient);
            w2 = gradient / (1 + gradient);
        end
        % Normalize weights to sum to 1
        sum_w = w1 + w2;
        w1 = w1 / sum_w;
        w2 = w2 / sum_w;
        if w2 < 1e-10
            ratio = "inf";
        else
            ratio = w1/w2;
        end
        % Print results
        fprintf('%4d | %10.4f | %10.6f | %9.4f | %12.4f | %10s\n', ...
            i, valid_cost(i), valid_volt(i), w1, w2, num2str(ratio));
    end
    [~, balanced_idx] = min(abs(arrayfun(@(i) (norm_cost(i) - norm_volt(i)), 1:length(valid_cost))));
    normalized_sum = norm_cost + norm_volt;
    [~, min_sum_idx] = min(normalized_sum);
    fprintf('\nRecommended Weight Sets:\n');
    fprintf('1. Most balanced weights: w1 = %.4f, w2 = %.4f (point %d)\n', ...
        1/(1+norm_volt(balanced_idx)/norm_cost(balanced_idx)), ...
        1/(1+norm_cost(balanced_idx)/norm_volt(balanced_idx)), balanced_idx);
     mid_idx = ceil(length(valid_cost) / 2);
    if mid_idx == balanced_idx
        fprintf('2. Middle Pareto point: w1 = %.4f, w2 = %.4f (point %d) (same as balanced)\n', ...
            1/(1+norm_volt(mid_idx)/norm_cost(mid_idx)), ...
            1/(1+norm_cost(mid_idx)/norm_volt(mid_idx)), mid_idx);
    else
        fprintf('2. Middle Pareto point: w1 = %.4f, w2 = %.4f (point %d)\n', ...
            1/(1+norm_volt(mid_idx)/norm_cost(mid_idx)), ...
            1/(1+norm_cost(mid_idx)/norm_volt(mid_idx)), mid_idx);
    end
    if min_sum_idx == balanced_idx || min_sum_idx == mid_idx
       else
        fprintf('3. Minimum normalized objective sum: w1 = %.4f, w2 = %.4f (point %d)\n', ...
            1/(1+norm_volt(min_sum_idx)/norm_cost(min_sum_idx)), ...
            1/(1+norm_cost(min_sum_idx)/norm_volt(min_sum_idx)), min_sum_idx);
    end
    figure('Name', 'Weight Distribution Along Pareto Front', 'Color', [1 1 1]);
    w1_values = zeros(length(valid_cost), 1);
    w2_values = zeros(length(valid_cost), 1);
    for i = 1:length(valid_cost)
        if i == 1
            dc = norm_cost(i+1) - norm_cost(i);
            dv = norm_volt(i+1) - norm_volt(i);
        elseif i == length(valid_cost)
            dc = norm_cost(i) - norm_cost(i-1);
            dv = norm_volt(i) - norm_volt(i-1);
        else
            dc = norm_cost(i+1) - norm_cost(i-1);
            dv = norm_volt(i+1) - norm_volt(i-1);
        end
       % Calculate weights
        if abs(dc) < 1e-10
            w1_values(i) = 0;
            w2_values(i) = 1;
        elseif abs(dv) < 1e-10
            w1_values(i) = 1;
            w2_values(i) = 0;
        else
            gradient = -dv/dc;
            w1_values(i) = 1 / (1 + gradient);
            w2_values(i) = gradient / (1 + gradient);
            % Normalize
            sum_w = w1_values(i) + w2_values(i);
            w1_values(i) = w1_values(i) / sum_w;
            w2_values(i) = w2_values(i) / sum_w;
        end
    end

    % Plotting weights
    plot(1:length(valid_cost), w1_values, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    hold on;
    plot(1:length(valid_cost), w2_values, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    plot([balanced_idx balanced_idx], [0 1], 'k--', 'LineWidth', 1.5);
    plot([mid_idx mid_idx], [0 1], 'g--', 'LineWidth', 1.5);
    grid on; box on;
    xlabel('Point Index on Pareto Front');
    ylabel('Weight Value');
    title('Weight Distribution Along Pareto Front');
    legend('w1 (Cost)', 'w2 (Voltage)', 'Balanced Point', 'Middle Point', 'Location', 'best');
    ylim([0 1]);
end
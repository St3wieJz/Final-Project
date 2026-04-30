%% Experiment 1: Intervention Timing vs. No Intervention
clear; clc;

% 1. Parameters (Fixed based on SIR-CIM framework)
p.alpha = 0.40; p.beta = 0.10; p.delta = 0.40; 
p.eta = 0.30;   p.mu = 0.05;   p.nu = 0.02; p.theta = 0.08;
N = 5000; T_max = 150; dt = 0.1;
delays = [5, 15, 30]; % Early, Mid, Late intervention
colors = {'#0072BD', '#D95319', '#EDB120'}; % Intervention colors

%% FIGURE 1: Spreader Density Comparison (S(t))
figure('Color', 'w', 'Name', 'Spreader Density Comparison');
hold on;

% --- Scenario 0: No Intervention Baseline ---
y0_baseline = [(N-1)/N; 1/N; 0; 0];
[t_base, y_base] = ode45(@(t,y) sir_cim_ode(t, y, p), 0:dt:T_max, y0_baseline);
plot(t_base, y_base(:,2), '--k', 'LineWidth', 2, 'DisplayName', 'No Intervention');

% --- Scenarios 1-3: With Intervention ---
for i = 1:length(delays)
    t_delay = delays(i);
    % Phase 1: Pre-intervention
    [t1, y1] = ode45(@(t,y) sir_cim_ode(t, y, p), 0:dt:t_delay, y0_baseline);
    % Phase 2: Post-intervention (Seed = 15/N)
    last_s = y1(end, :);
    y0_p2 = [max(0, last_s(1)-15/N), last_s(2), 15/N, last_s(4)];
    [t2, y2] = ode45(@(t,y) sir_cim_ode(t, y, p), t_delay:dt:T_max, y0_p2);
    
    plot([t1; t2], [y1(:,2); y2(:,2)], 'Color', colors{i}, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Intervention at t=%d', t_delay));
end

title('Spreader Density S(t): Impact of Truth-Spreader Timing', 'FontSize', 14);
xlabel('Time (t)'); ylabel('Density of Spreaders S(t)');
grid on; legend('Location', 'northeast');
set(gca, 'FontSize', 11);

%% FIGURE 2: Full System Dynamics (C, S, T, R) Comparison
% This figure directly visualizes the "Rumor Size" (Final R)
fig2 = figure('Color', 'w', 'Position', [100 100 1200 800]);
tiledlayout(2, 2, 'TileSpacing', 'Compact');

% We will plot: [No Intervention, t=5, t=15, t=30]
all_delays = [Inf, 5, 15, 30]; 
titles = {'No Intervention', 'Early Intervention (t=5)', ...
          'Mid Intervention (t=15)', 'Late Intervention (t=30)'};

for k = 1:length(all_delays)
    nexttile;
    if isinf(all_delays(k))
        t_plot = t_base; y_plot = y_base;
    else
        t_d = all_delays(k);
        [t1, y1] = ode45(@(t,y) sir_cim_ode(t, y, p), 0:dt:t_d, y0_baseline);
        y0_p2 = [max(0, y1(end,1)-15/N), y1(end,2), 15/N, y1(end,4)];
        [t2, y2] = ode45(@(t,y) sir_cim_ode(t, y, p), t_d:dt:T_max, y0_p2);
        t_plot = [t1; t2]; y_plot = [y1; y2];
        xline(t_d, ':k', 'HandleVisibility', 'off'); % Marker for intervention
    end
    
    plot(t_plot, y_plot(:,1), 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); hold on; % C
    plot(t_plot, y_plot(:,2), 'r', 'LineWidth', 1.5);                             % S
    plot(t_plot, y_plot(:,3), 'b', 'LineWidth', 1.5);                             % T
    plot(t_plot, y_plot(:,4), 'g', 'LineWidth', 1.5);                             % R
    
    title(titles{k});
    if k >= 3, xlabel('Time'); end
    if mod(k,2) ~= 0, ylabel('Density'); end
    grid on; ylim([0 1]);
    if k == 1, legend({'C', 'S', 'T', 'R'}, 'Location', 'northeast'); end
end

% --- ODE Function ---
function dydt = sir_cim_ode(t, y, p)
    C = y(1); S = y(2); T = y(3); R = y(4);
    % dC/dt = -(alpha + beta)CS - delta*CT
    dCdt = -(p.alpha + p.beta) * C * S - p.delta * C * T;
    % dSdt = alpha*CS - mu*S - theta*S(S + R) - eta*ST
    dSdt = p.alpha * C * S - p.mu * S - p.theta * S * (S + R) - p.eta * S * T;
    % dTdt = delta*CT - nu*T
    dTdt = p.delta * C * T - p.nu * T;
    % dRdt = beta*CS + mu*S + theta*S(S + R) + nu*T + eta*ST
    dRdt = p.beta * C * S + p.mu * S + p.theta * S * (S + R) + p.nu * T + p.eta * S * T;
    dydt = [dCdt; dSdt; dTdt; dRdt];
end
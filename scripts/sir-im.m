
% Reproduce the core figures for the SIR-IM rumor spreading model.

clear; clc; close all;

%% 1. Global Parameters and Initial Conditions
N = 5000;               % Size of the network 
tspan = [0 100];        % Simulation time span
y0 = [(N-1)/N; 1/N; 0]; % Initial densities: [Credulous; Spreader; Stifler]

%% 2. Recreating Figure 4d: Baseline SIR-IM Model
% Default parameters for baseline
p_base = 0.2; 
theta_base = 0.02; 
lambda_base = 0.1; 
beta2_base = 0.2;

[t, y] = ode45(@(t,y) sir_im_ode(t, y, p_base, theta_base, lambda_base, beta2_base, N), tspan, y0);

figure('Name', 'Figure 4d: SIR-IM Baseline');
plot(t, y(:,1), 'r-+', 'LineWidth', 1.5); hold on;
plot(t, y(:,2), 'g-', 'LineWidth', 1.5);
plot(t, y(:,3), 'k-*', 'LineWidth', 1.5);
title('SIR-IM Model (Baseline)');
xlabel('Time'); ylabel('Density');
legend('Credulous', 'Spreader', 'Stifler');
grid on;

%% 3. Recreating Figure 5: Effect of individual opinion (P)
P_values = [0.1, 0.3, 0.9];
figure('Name', 'Figure 5: Effect of P', 'Position', [100 100 1200 400]);

for i = 1:length(P_values)
    [t_P, y_P] = ode45(@(t,y) sir_im_ode(t, y, P_values(i), theta_base, lambda_base, beta2_base, N), tspan, y0);
    
    subplot(1, 3, i);
    plot(t_P, y_P(:,1), 'r-+', 'LineWidth', 1.5); hold on;
    plot(t_P, y_P(:,2), 'g-', 'LineWidth', 1.5);
    plot(t_P, y_P(:,3), 'k-*', 'LineWidth', 1.5);
    title(sprintf('P = %.1f', P_values(i)));
    xlabel('Time'); ylabel('Density');
    legend('Credulous', 'Spreader', 'Stifler');
    grid on;
end

%% 4. Recreating Figure 6: Effect of forgetting rate parameter (\lambda)
lambda_values = [0.1, 0.3, 0.9];
figure('Name', 'Figure 6: Effect of \lambda');
colors = {'r-o', 'g-', 'k-+'};

for i = 1:length(lambda_values)
    [t_L, y_L] = ode45(@(t,y) sir_im_ode(t, y, p_base, theta_base, lambda_values(i), beta2_base, N), tspan, y0);
    plot(t_L, y_L(:,2), colors{i}, 'LineWidth', 1.5); hold on;
end
title('The density of spreader with different \lambda');
xlabel('Time'); ylabel('Density of Spreader');
legend('\lambda = 0.1', '\lambda = 0.3', '\lambda = 0.9');
grid on;

%% 5. Recreating Figure 7: Effect of immune probability (\theta)
theta_values = [0.02, 0.1, 0.5];
figure('Name', 'Figure 7: Effect of \theta');
colors = {'r-o', 'g-', 'k-+'};

for i = 1:length(theta_values)
    [t_T, y_T] = ode45(@(t,y) sir_im_ode(t, y, p_base, theta_values(i), lambda_base, beta2_base, N), tspan, y0);
    plot(t_T, y_T(:,2), colors{i}, 'LineWidth', 1.5); hold on;
end
title('The density of spreader with different \theta');
xlabel('Time'); ylabel('Density of Spreader');
legend('\theta = 0.02', '\theta = 0.1', '\theta = 0.5');
grid on;

%% 6. Recreating Figure 8: Effect of stifling rate (\beta_2)
beta2_values = [0.2, 0.4, 0.8];
figure('Name', 'Figure 8: Effect of \beta_2', 'Position', [100 100 1200 400]);
colors = {'r-o', 'g-', 'k-+'};

% Compute all three runs first
data_beta = cell(3,2);
for i = 1:length(beta2_values)
    [t_B, y_B] = ode45(@(t,y) sir_im_ode(t, y, p_base, theta_base, lambda_base, beta2_values(i), N), tspan, y0);
    data_beta{i, 1} = t_B;
    data_beta{i, 2} = y_B;
end

% Plotting the 3 subplots (Credulous, Spreader, Stifler)
titles = {'The credulous with different \beta_2', 'The spreader with different \beta_2', 'The stifler with different \beta_2'};
y_labels = {'The density of credulous', 'The density of spreader', 'The density of stifler'};

for state_idx = 1:3 % 1=C, 2=S, 3=R
    subplot(1, 3, state_idx);
    for i = 1:length(beta2_values)
        plot(data_beta{i,1}, data_beta{i,2}(:, state_idx), colors{i}, 'LineWidth', 1.5); hold on;
    end
    title(titles{state_idx});
    xlabel('Time'); ylabel(y_labels{state_idx});
    legend('\beta_2 = 0.2', '\beta_2 = 0.4', '\beta_2 = 0.8');
    grid on;
end


function dydt = sir_im_ode(t, y, p, theta, lambda, beta2, N)
    % Extract current densities
    C = y(1);
    S = y(2);
    R = y(3);
    
    % Calculate the current number of spreaders (m)
    m = S * N; 
    
    % Calculate dynamic spreading and skeptical probabilities
    alpha1_prime = 1 - (1 - p)^m;
    alpha2_prime = (1 - p)^m;
    
    alpha1 = alpha1_prime / (alpha1_prime + alpha2_prime + p);
    alpha2 = alpha2_prime / (alpha1_prime + alpha2_prime + p);
    
    % Calculate dynamic forgetting rate
    beta1 = 1 - exp(-lambda * t);
    
    % --- THE MISSING LINK ---
    % The paper generated figures using a Barabasi-Albert (BA) network simulation, 
    % which has an average degree (k) > 1. The ODEs in the paper assumed a contact rate of 1.
    % We add an effective network multiplier (k_avg) to bridge this gap.
    k_avg = 3; 
    
    % Apply the network multiplier to all contact-based interactions
    dCdt = -k_avg * (alpha1 + alpha2 + theta) * C * S;
    dSdt = k_avg * alpha1 * C * S - beta1 * S - k_avg * beta2 * (S + R) * S;
    dRdt = k_avg * (alpha2 + theta) * C * S + beta1 * S + k_avg * beta2 * (S + R) * S;
    
    % Return column vector
    dydt = [dCdt; dSdt; dRdt];
end
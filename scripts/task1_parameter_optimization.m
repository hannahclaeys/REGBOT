%% REGBOT Wheel Speed Controller - Parameter Optimization
% This script loads real REGBOT data, identifies the transfer function,
% and optimizes PI controller parameters by varying gamma_M and N_i.
%
% The script tests different combinations to find which configuration gives
% the most stable response (lowest overshoot and fastest settling time).

clear all; close all; clc;

%% =========================================================
% 1. LOAD AND PROCESS DATA (from REGBOT_script.m)
% ==========================================================

fprintf('Loading REGBOT measurement data...\n');

data = readtable('W6_BACK_2sec_code_Kim.txt');
data = fillmissing(data, 'nearest');

% Remove false start
data = data(10:end, :);

% Extract signals
t = table2array(data(:, 1));
u_L = table2array(data(:, 8));  % Motor voltage left
u_R = table2array(data(:, 9));  % Motor voltage right
v_L = table2array(data(:, 10)); % Wheel velocity left
v_R = table2array(data(:, 11)); % Wheel velocity right

% Sampling time
T_s = t(2) - t(1);

% Trim data
idx = max(find(u_R > 3, 1, 'first') - 2000, 1);
t = t(idx:end);
u_L = u_L(idx:end);
u_R = u_R(idx:end);
v_L = v_L(idx:end);
v_R = v_R(idx:end);

% Start time at 0
t = t - t(1);

% Remove offsets from velocity
v_L = v_L - v_L(1);
v_R = v_R - v_R(1);

% Remove outliers
bad_L = (t > 1) & (v_L < 0.1);
bad_R = (t > 1) & (v_R < 0.1);
v_L(bad_L) = NaN;
v_R(bad_R) = NaN;
v_L = fillmissing(v_L, 'linear');
v_R = fillmissing(v_R, 'linear');

% Average signals
volt = 0.5 * (u_L + u_R);
vel = 0.5 * (v_L + v_R);

fprintf('✓ Data loaded successfully (%.3f seconds of measurement)\n', t(end));

%% =========================================================
% 2. IDENTIFY TRANSFER FUNCTION
% ==========================================================

fprintf('\nIdentifying transfer function G_avg(s) from measurement data...\n');

idd_avg = iddata(vel, volt, T_s);
G_avg = tfest(idd_avg, 1, 0);  % 1 pole, 0 zeros

[num_G, den_G] = tfdata(G_avg, 'v');
fprintf('✓ Transfer function identified:\n');
fprintf('  G_avg(s) = %.4f / (s + %.2f)\n', num_G(1), -den_G(2));

% Extract plant parameters
K_plant = num_G(1);
pole_plant = -den_G(2);

%% =========================================================
% 3. DEFINE PARAMETER RANGES FOR OPTIMIZATION
% ==========================================================

fprintf('\n=========================================================\n');
fprintf('PARAMETER OPTIMIZATION\n');
fprintf('=========================================================\n');

% Design parameters to optimize
gamma_M_range = [45, 60, 75];   % Phase margin in degrees
N_i_range = [2, 3, 4, 5];       % Filter ratio

% Reference velocity
ref_velocity = 0.5;  % m/s

%% =========================================================
% 4. LOOP THROUGH ALL COMBINATIONS
% ==========================================================

results = struct();
idx = 1;

% Create figure for step responses
fig_steps = figure('Position', [100 100 1400 900]);
fig_steps.Name = 'PI Controller Parameter Optimization - Step Responses';

for gamma_M = gamma_M_range
    for N_i = N_i_range
        
        % 1. Compute phase contribution from integrator
        phi_i_rad = -atan(1/N_i);
        phi_i_deg = rad2deg(phi_i_rad);
        
        % 2. Find crossover frequency omega_c from Bode plot
        % Phase equation: angle(G(j*omega_c)) = gamma_M - phi_i - 180
        target_phase = gamma_M + phi_i_deg - 180;  % in degrees
        
        % Create frequency vector and find where phase matches
        w = logspace(-1, 2, 10000);
        [mag, phase] = bode(G_avg, w);
        mag = squeeze(mag);
        phase = squeeze(phase);
        
        % Find closest frequency to target phase
        [~, idx_w] = min(abs(phase - target_phase));
        omega_c = w(idx_w);
        
        % 3. Compute integral time constant
        tau_i = N_i / omega_c;
        
        % 4. Compute proportional gain K_P
        % |C_PI(j*omega_c) * G(j*omega_c)| = 1
        s_c = 1j * omega_c;
        C_PI_mag = abs((tau_i * s_c + 1) / (tau_i * s_c));
        G_mag = abs(K_plant / (s_c + pole_plant));
        K_P = 1 / (C_PI_mag * G_mag);
        
        % 5. Create PI controller transfer function
        C_PI = tf([K_P * tau_i, K_P], [tau_i, 0]);
        
        % 6. Compute closed-loop system
        sys_openloop = C_PI * G_avg;
        sys_closedloop = feedback(sys_openloop, 1);
        
        % Store results
        results(idx).gamma_M = gamma_M;
        results(idx).N_i = N_i;
        results(idx).phi_i_deg = phi_i_deg;
        results(idx).omega_c = omega_c;
        results(idx).tau_i = tau_i;
        results(idx).K_P = K_P;
        results(idx).C_PI = C_PI;
        results(idx).sys_closedloop = sys_closedloop;
        
        % Print parameters
        fprintf('\n--- Configuration %d ---\n', idx);
        fprintf('γ_M = %6.2f°   |   N_i = %d\n', gamma_M, N_i);
        fprintf('φ_i = %7.2f°   |   ω_c = %6.3f rad/s\n', phi_i_deg, omega_c);
        fprintf('τ_i = %7.4f   |   K_P = %7.4f\n', tau_i, K_P);
        
        idx = idx + 1;
    end
end

%% =========================================================
% 5. COMPUTE PERFORMANCE METRICS (Overshoot & Settling Time)
% ==========================================================

fprintf('\n\nComputing performance metrics...\n');

t_sim = 0:0.01:5;
performance = struct();

for i = 1:length(results)
    % Get step response
    y_sim = step(results(i).sys_closedloop, t_sim);
    
    % Compute overshoot
    max_response = max(y_sim);
    overshoot_percent = (max_response - ref_velocity) / ref_velocity * 100;
    
    % Compute settling time (2% criterion manually)
    steady_state = y_sim(end);
    tolerance_band = 0.02 * steady_state;
    
    settled_idx = find(abs(y_sim - steady_state) > tolerance_band, 1, 'last');
    if isempty(settled_idx)
        settling_time = t_sim(1);
    else
        settling_time = t_sim(settled_idx);
    end
    
    % Store metrics
    performance(i).overshoot = overshoot_percent;
    performance(i).settling_time = settling_time;
    performance(i).max_response = max_response;
end

%% =========================================================
% 6. PLOT STEP RESPONSES FOR ALL CONFIGURATIONS
% ==========================================================

fprintf('Generating step response plots...\n\n');

figure(fig_steps);
num_configs = length(results);
rows = ceil(sqrt(num_configs));
cols = ceil(num_configs / rows);

for i = 1:num_configs
    subplot(rows, cols, i);
    
    % Step response
    t_sim = 0:0.01:5;
    y_sim = step(results(i).sys_closedloop, t_sim);
    
    plot(t_sim, y_sim, 'b-', 'LineWidth', 2.5);
    hold on;
    plot(t_sim, ones(size(t_sim)) * ref_velocity, 'r--', 'LineWidth', 1.5);
    grid on;
    
    title(sprintf('γ_M = %.0f°, N_i = %d\nK_P = %.3f, τ_i = %.4f', ...
                  results(i).gamma_M, results(i).N_i, ...
                  results(i).K_P, results(i).tau_i), 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    ylim([0 0.7]);
    set(gca, 'FontSize', 9);
    
    % Get performance metrics
    overshoot_percent = performance(i).overshoot;
    settling_time = performance(i).settling_time;
    
    % Color background based on performance
    if overshoot_percent < 10
        bg_color = [0.9 1.0 0.9];  % Light green (good)
    elseif overshoot_percent < 20
        bg_color = [1.0 1.0 0.9];  % Light yellow (ok)
    else
        bg_color = [1.0 0.9 0.9];  % Light red (poor)
    end
    set(gca, 'Color', bg_color);
    
    text(0.5, 0.05, sprintf('OS = %.1f%%\nT_s = %.2fs', overshoot_percent, settling_time), ...
         'Units', 'normalized', 'FontSize', 9, 'FontWeight', 'bold',...
         'BackgroundColor', 'white', 'EdgeColor', 'black');
end

sgtitle('PI Controller Parameter Optimization - Step Response Comparison', ...
        'FontSize', 14, 'FontWeight', 'bold');

%% =========================================================
% 7. SUMMARY TABLE
% ==========================================================

fprintf('==================================================== SUMMARY TABLE ====================================================\n');
fprintf('%5s | %8s | %5s | %10s | %10s | %10s | %8s | %8s\n', ...
        'Idx', 'gamma_M', 'N_i', 'omega_c', 'tau_i', 'K_P', 'Oversh.', 'Settling');
fprintf('-----------------------------------------------------------------------------------------------------------------------\n');

for i = 1:num_configs
    fprintf('%5d | %8.1f° | %5d | %10.3f | %10.4f | %10.4f | %7.2f%% | %8.2f\n', ...
            i, results(i).gamma_M, results(i).N_i, results(i).omega_c, ...
            results(i).tau_i, results(i).K_P, performance(i).overshoot, performance(i).settling_time);
end

fprintf('-----------------------------------------------------------------------------------------------------------------------\n');

%% =========================================================
% 8. FIND BEST CONFIGURATION
% ==========================================================

fprintf('\nSearching for optimal configuration based on stability criteria...\n');
best_score = inf;
best_idx = 1;

for i = 1:num_configs
    % Score: minimize overshoot (heavily weighted) and settling time
    % Overshoot is more important for stability
    score = 2 * performance(i).overshoot + 0.5 * performance(i).settling_time;
    
    if score < best_score
        best_score = score;
        best_idx = i;
    end
end

fprintf('\n');
fprintf('╔══════════════════════════════════════════════╗\n');
fprintf('║        ⭐ OPTIMAL CONFIGURATION FOUND ⭐       ║\n');
fprintf('╚══════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('  Configuration Index: %d\n', best_idx);
fprintf('  Phase Margin (γ_M): %.1f°\n', results(best_idx).gamma_M);
fprintf('  Filter Ratio (N_i):  %d\n', results(best_idx).N_i);
fprintf('  Proportional Gain (K_P): %.4f\n', results(best_idx).K_P);
fprintf('  Integral Time (τ_i):     %.4f\n', results(best_idx).tau_i);
fprintf('  Crossover Freq (ω_c):    %.4f rad/s\n', results(best_idx).omega_c);
fprintf('\n');
fprintf('  Performance:\n');
fprintf('    • Overshoot: %.2f%%\n', performance(best_idx).overshoot);
fprintf('    • Settling Time: %.2f s\n', performance(best_idx).settling_time);
fprintf('\n');

%% =========================================================
% 9. PLOT BODE DIAGRAMS FOR BEST CONFIGURATION
% ==========================================================

fprintf('Generating Bode plots for best configuration...\n');

fig_bode = figure('Position', [100 100 1200 600]);
fig_bode.Name = sprintf('Best Configuration (gamma_M=%.0f°, N_i=%d) - Bode Plot', ...
                        results(best_idx).gamma_M, results(best_idx).N_i);

% Open-loop Bode
subplot(1, 2, 1);
[mag, phase, w] = bode(results(best_idx).C_PI * G_avg);
semilogx(w, squeeze(phase), 'b-', 'LineWidth', 2.5);
grid on;
xlabel('Frequency (rad/s)', 'FontSize', 11);
ylabel('Phase (degrees)', 'FontSize', 11);
title('Open-Loop Bode Plot - Phase', 'FontSize', 12, 'FontWeight', 'bold');
hold on;
yline(-180, 'r--', 'LineWidth', 1.5, 'DisplayName', '-180°');
legend('Phase response', 'Location', 'best');
set(gca, 'FontSize', 10);

% Magnitude
subplot(1, 2, 2);
semilogx(w, 20*log10(squeeze(mag)), 'b-', 'LineWidth', 2.5);
grid on;
xlabel('Frequency (rad/s)', 'FontSize', 11);
ylabel('Magnitude (dB)', 'FontSize', 11);
title('Open-Loop Bode Plot - Magnitude', 'FontSize', 12, 'FontWeight', 'bold');
hold on;
yline(0, 'r--', 'LineWidth', 1.5, 'DisplayName', '0 dB');
legend('Magnitude response', 'Location', 'best');
set(gca, 'FontSize', 10);

sgtitle(sprintf('Best Configuration Analysis (γ_M=%.0f°, N_i=%d)', ...
                results(best_idx).gamma_M, results(best_idx).N_i), ...
        'FontSize', 14, 'FontWeight', 'bold');

%% =========================================================
% 10. COMPARE WITH MEASURED DATA
% ==========================================================

fprintf('Generating comparison with measured data...\n');

% Step response of best controller on identified model
[y_model, t_model] = step(results(best_idx).sys_closedloop, t(end));

fig_compare = figure('Position', [100 100 1200 600]);
fig_compare.Name = 'Model vs Measured Response';

% Plot comparison
plot(t_model, y_model, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Model (Linear)');
hold on;
plot(t, vel, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Measured Data');
grid on;

xlabel('Time (s)', 'FontSize', 11);
ylabel('Velocity (m/s)', 'FontSize', 11);
title(sprintf('Step Response Comparison\n(γ_M=%.0f°, N_i=%d, K_P=%.4f)', ...
             results(best_idx).gamma_M, results(best_idx).N_i, results(best_idx).K_P), ...
      'FontSize', 12, 'FontWeight', 'bold');
legend('Model Response', 'Measured Response', 'Location', 'best', 'FontSize', 11);
set(gca, 'FontSize', 10);

fprintf('✓ Comparison plot generated\n');
fprintf('\nNote: The linear model does not capture nonlinear effects like tilt coupling.\n');
fprintf('For real robot testing, implement the controller and validate on REGBOT hardware.\n');

%% =========================================================
% 11. SAVE OPTIMAL PARAMETERS
% ==========================================================

fprintf('\n=========================================================\n');
fprintf('SAVING RESULTS\n');
fprintf('=========================================================\n');

optimal_params.gamma_M = results(best_idx).gamma_M;
optimal_params.N_i = results(best_idx).N_i;
optimal_params.K_P = results(best_idx).K_P;
optimal_params.tau_i = results(best_idx).tau_i;
optimal_params.omega_c = results(best_idx).omega_c;
optimal_params.phi_i = results(best_idx).phi_i_deg;
optimal_params.plant_K = K_plant;
optimal_params.plant_pole = pole_plant;
optimal_params.performance_overshoot = performance(best_idx).overshoot;
optimal_params.performance_settling_time = performance(best_idx).settling_time;

save('optimal_PI_parameters.mat', 'optimal_params', 'results', 'performance');
fprintf('✓ Results saved to: optimal_PI_parameters.mat\n');

fprintf('\n✅ OPTIMIZATION COMPLETE!\n\n');

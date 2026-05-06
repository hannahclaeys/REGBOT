%% REGBOT Wheel Velocity Controller - Simulink Parameter Optimization
% This script optimizes PI controller parameters (gamma_M, N_i) by running
% the nonlinear Simulink model (regbot_1mg_task2.slx) and measuring stability
% metrics (overshoot and settling time) on wheel velocity (lin_vel).
%
% The script tests a grid of parameters and finds the optimal configuration.

clear all; close all; clc;

%% =========================================================
% 1. LOAD MEASURED DATA & IDENTIFY TRANSFER FUNCTION
% ==========================================================

fprintf('╔════════════════════════════════════════════════════════╗\n');
fprintf('║  REGBOT Wheel Velocity Controller - Optimization      ║\n');
fprintf('╚════════════════════════════════════════════════════════╝\n\n');

fprintf('Loading REGBOT measurement data...\n');

data = readtable('W6_BACK_2sec_code_Kim.txt');
data = fillmissing(data, 'nearest');
data = data(10:end, :);

t = table2array(data(:, 1));
u_L = table2array(data(:, 8));
u_R = table2array(data(:, 9));
v_L = table2array(data(:, 10));
v_R = table2array(data(:, 11));

T_s = t(2) - t(1);
idx = max(find(u_R > 3, 1, 'first') - 2000, 1);

t = t(idx:end);
u_L = u_L(idx:end);
u_R = u_R(idx:end);
v_L = v_L(idx:end);
v_R = v_R(idx:end);

t = t - t(1);
v_L = v_L - v_L(1);
v_R = v_R - v_R(1);

bad_L = (t > 1) & (v_L < 0.1);
bad_R = (t > 1) & (v_R < 0.1);
v_L(bad_L) = NaN;
v_R(bad_R) = NaN;
v_L = fillmissing(v_L, 'linear');
v_R = fillmissing(v_R, 'linear');

volt = 0.5 * (u_L + u_R);
vel = 0.5 * (v_L + v_R);

fprintf('✓ Data loaded successfully\n');

%% =========================================================
% 2. IDENTIFY TRANSFER FUNCTION
% ==========================================================

fprintf('Identifying transfer function G_avg(s)...\n');

idd_avg = iddata(vel, volt, T_s);
G_avg = tfest(idd_avg, 1, 0);

[num_G, den_G] = tfdata(G_avg, 'v');
K_plant = num_G(1);
pole_plant = -den_G(2);

fprintf('✓ G_avg(s) = %.4f / (s + %.2f)\n\n', K_plant, pole_plant);

%% =========================================================
% 3. DEFINE PARAMETER RANGES
% ==========================================================

gamma_M_range = 40:1:70;        % 40° to 70° in steps of 1°
N_i_range = 1:0.5:5;            % 1 to 5 in steps of 0.5

fprintf('Parameter ranges:\n');
fprintf('  gamma_M: [%.0f:1:%.0f] → %d values\n', gamma_M_range(1), gamma_M_range(end), length(gamma_M_range));
fprintf('  N_i:     [%.1f:0.5:%.1f] → %d values\n', N_i_range(1), N_i_range(end), length(N_i_range));
fprintf('  Total configurations: %d\n\n', length(gamma_M_range) * length(N_i_range));

%% =========================================================
% 4. CALCULATE CONTROLLER PARAMETERS FOR ALL CONFIGURATIONS
% ==========================================================

fprintf('Calculating controller parameters for all configurations...\n');

idx = 1;
for gamma_M = gamma_M_range
    for N_i = N_i_range
        
        % Phase contribution from integrator
        phi_i_rad = -atan(1/N_i);
        phi_i_deg = rad2deg(phi_i_rad);
        
        % Target phase for crossover frequency
        target_phase = gamma_M + phi_i_deg - 180;
        
        % Find omega_c from Bode plot
        w = logspace(-1, 2, 10000);
        [mag, phase] = bode(G_avg, w);
        mag = squeeze(mag);
        phase = squeeze(phase);
        
        [~, idx_w] = min(abs(phase - target_phase));
        omega_c = w(idx_w);
        
        % Integral time constant
        tau_i = N_i / omega_c;
        
        % Proportional gain K_P
        s_c = 1j * omega_c;
        C_PI_mag = abs((tau_i * s_c + 1) / (tau_i * s_c));
        G_mag = abs(K_plant / (s_c + pole_plant));
        K_P = 1 / (C_PI_mag * G_mag);
        
        % Store parameters
        params(idx).gamma_M = gamma_M;
        params(idx).N_i = N_i;
        params(idx).phi_i = phi_i_deg;
        params(idx).omega_c = omega_c;
        params(idx).tau_i = tau_i;
        params(idx).K_P = K_P;
        
        idx = idx + 1;
    end
end

fprintf('✓ %d configurations ready\n\n', idx-1);

%% =========================================================
% 5. LOAD SIMULINK MODEL & RUN SIMULATIONS
% ==========================================================

fprintf('═════════════════════════════════════════════════════════\n');
fprintf('RUNNING SIMULINK SIMULATIONS\n');
fprintf('═════════════════════════════════════════════════════════\n\n');

model_name = 'regbot_1mg_task2';
sim_time = 420;  % 7 minutes in seconds

% Load model
load_system(model_name);

fprintf('Simulating %d configurations (this may take a while)...\n\n', idx-1);

performance = struct();
num_configs = idx - 1;

% Progress bar setup
start_time = tic;

for i = 1:num_configs
    
    % Extract parameters
    gamma_M = params(i).gamma_M;
    N_i = params(i).N_i;
    K_P = params(i).K_P;
    tau_i = params(i).tau_i;
    
    % Set Simulink parameters
    set_param(model_name, 'StopTime', num2str(sim_time));
    
    % Set controller gains in model
    set_param([model_name '/Kpwv'], 'Gain', num2str(K_P));
    set_param([model_name '/1//tiwv'], 'Gain', num2str(1/tau_i));
    
    % Run simulation
    try
        sim_output = sim(model_name, 'StopTime', num2str(sim_time));
        
        % Extract wheel velocity from simulation output
        % The model outputs time series data - we need to find lin_vel
        if isfield(sim_output, 'lin_vel')
            lin_vel_data = sim_output.lin_vel;
            time_data = sim_output.tout;
        else
            % Try to get from workspace (Scope outputs to base)
            lin_vel_data = evalin('base', 'lin_vel_data.signals.values', []);
            time_data = evalin('base', 'lin_vel_data.time', []);
        end
        
        if isempty(lin_vel_data)
            error('Could not find lin_vel output from simulation');
        end
        
        % Ensure column vector
        if size(lin_vel_data, 2) > 1
            lin_vel_data = lin_vel_data(:, 1);
        end
        
        % Calculate performance metrics
        ref_velocity = 0.5;  % Reference velocity
        
        % Max overshoot (%)
        max_vel = max(lin_vel_data);
        overshoot = (max_vel - ref_velocity) / ref_velocity * 100;
        
        % Settling time (2% criterion)
        steady_state = lin_vel_data(end);
        tolerance_band = 0.02 * ref_velocity;
        
        settled_idx = find(abs(lin_vel_data - ref_velocity) > tolerance_band, 1, 'last');
        if isempty(settled_idx)
            settling_time = time_data(1);
        else
            settling_time = time_data(min(settled_idx + 1, length(time_data)));
        end
        
        % Store results
        performance(i).overshoot = overshoot;
        performance(i).settling_time = settling_time;
        performance(i).success = true;
        
    catch ME
        % Simulation failed
        performance(i).overshoot = inf;
        performance(i).settling_time = inf;
        performance(i).success = false;
        performance(i).error = ME.message;
    end
    
    % Progress output
    if mod(i, 5) == 0 || i == 1
        elapsed = toc(start_time);
        rate = i / elapsed;
        remaining = (num_configs - i) / rate;
        fprintf('[%3d/%3d] γ_M=%.0f° N_i=%.1f | OS=%6.2f%% T_s=%6.2f | ETA: %.1f min\n', ...
                i, num_configs, gamma_M, N_i, performance(i).overshoot, ...
                performance(i).settling_time, remaining/60);
    end
end

% Unload model
close_system(model_name, 0);

fprintf('\n✓ All simulations completed!\n\n');

%% =========================================================
% 6. ANALYZE RESULTS & DISPLAY TABLE
% ==========================================================

fprintf('╔════════════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                          RESULTS TABLE                                   ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════════════╝\n\n');

fprintf('%4s | %8s | %5s | %10s | %10s | %8s | %10s\n', ...
        'Idx', 'gamma_M', 'N_i', 'K_P', 'tau_i', 'Oversh.', 'Settling');
fprintf('────────────────────────────────────────────────────────────────────────────\n');

best_score = inf;
best_idx = 1;

for i = 1:num_configs
    if performance(i).success
        gamma_M = params(i).gamma_M;
        N_i = params(i).N_i;
        K_P = params(i).K_P;
        tau_i = params(i).tau_i;
        
        fprintf('%4d | %8.1f° | %5.1f | %10.4f | %10.4f | %7.2f%% | %9.2f s\n', ...
                i, gamma_M, N_i, K_P, tau_i, performance(i).overshoot, performance(i).settling_time);
        
        % Score calculation: minimize overshoot (weight 3) and settling time (weight 1)
        score = 3 * performance(i).overshoot + performance(i).settling_time;
        
        if score < best_score
            best_score = score;
            best_idx = i;
        end
    else
        fprintf('%4d | (simulation failed)\n', i);
    end
end

fprintf('────────────────────────────────────────────────────────────────────────────\n\n');

%% =========================================================
% 7. DISPLAY OPTIMAL CONFIGURATION
% ==========================================================

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════╗\n');
fprintf('║           ⭐ OPTIMAL CONFIGURATION FOUND ⭐             ║\n');
fprintf('╚════════════════════════════════════════════════════════╝\n\n');

fprintf('  Configuration Index:        %d\n', best_idx);
fprintf('  Phase Margin (γ_M):         %.1f°\n', params(best_idx).gamma_M);
fprintf('  Filter Ratio (N_i):         %.1f\n', params(best_idx).N_i);
fprintf('\n');
fprintf('  Controller Parameters:\n');
fprintf('    • Proportional Gain (K_P): %.4f\n', params(best_idx).K_P);
fprintf('    • Integral Time (τ_i):     %.4f\n', params(best_idx).tau_i);
fprintf('    • Phase Contribution (φ_i): %.2f°\n', params(best_idx).phi_i);
fprintf('    • Crossover Freq (ω_c):    %.4f rad/s\n', params(best_idx).omega_c);
fprintf('\n');
fprintf('  Performance Metrics:\n');
fprintf('    • Max Overshoot:           %.2f%%\n', performance(best_idx).overshoot);
fprintf('    • Settling Time (2%):      %.2f s\n', performance(best_idx).settling_time);
fprintf('\n');

%% =========================================================
% 8. SAVE RESULTS
% ==========================================================

optimal_result.gamma_M = params(best_idx).gamma_M;
optimal_result.N_i = params(best_idx).N_i;
optimal_result.K_P = params(best_idx).K_P;
optimal_result.tau_i = params(best_idx).tau_i;
optimal_result.phi_i = params(best_idx).phi_i;
optimal_result.omega_c = params(best_idx).omega_c;
optimal_result.overshoot = performance(best_idx).overshoot;
optimal_result.settling_time = performance(best_idx).settling_time;

save('optimal_PI_controller_simulink.mat', 'optimal_result', 'params', 'performance');
fprintf('✓ Results saved to: optimal_PI_controller_simulink.mat\n\n');

fprintf('╔════════════════════════════════════════════════════════╗\n');
fprintf('║                ✅ OPTIMIZATION COMPLETE              ║\n');
fprintf('╚════════════════════════════════════════════════════════╝\n\n');

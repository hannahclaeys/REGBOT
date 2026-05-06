%% =========================================================
% REGBOT
% ==========================================================

clear
close all
clc

%% =========================================================
% 1. LOAD DATA
% ==========================================================

data = readtable('W6_BACK_2sec_code_Kim.txt');
data = fillmissing(data,'nearest');

% Remove false start
data = data(10:end,:);

% Time
t = table2array(data(:,1));

% Input signals (motor voltage)
u_L = table2array(data(:,8));
u_R = table2array(data(:,9));

% Output signals (wheel velocity)
v_L = table2array(data(:,10));
v_R = table2array(data(:,11));

% Sampling time
T_s = t(2) - t(1);

%% =========================================================
% 2. TRIM DATA
% Keep data from just before the voltage changes from 3V to 4V
% ==========================================================

idx = max(find(u_R > 3,1,'first') - 2000, 1);

t   = t(idx:end);
u_L = u_L(idx:end);
u_R = u_R(idx:end);
v_L = v_L(idx:end);
v_R = v_R(idx:end);

% Start time at 0
t = t - t(1);

%% =========================================================
% 3. REMOVE OFFSETS
% Needed if we want to use step(...)
% ==========================================================

%u_L = u_L - u_L(1);
%u_R = u_R - u_R(1);

v_L = v_L - v_L(1);
v_R = v_R - v_R(1);

%% =========================================================
% 4. REMOVE OUTLIERS MANUALLY
% Remove clearly wrong low spikes
% ==========================================================

bad_L = (t > 1) & (v_L < 0.1);
bad_R = (t > 1) & (v_R < 0.1);

v_L(bad_L) = NaN;
v_R(bad_R) = NaN;

v_L = fillmissing(v_L,'linear');
v_R = fillmissing(v_R,'linear');

%% =========================================================
% 5. MAKE AVERAGE SIGNALS
% ==========================================================

volt = 0.5*(u_L + u_R);
vel  = 0.5*(v_L + v_R);

%% =========================================================
% 6. PLOT THE MEASURED DATA
% ==========================================================

figure(100)

subplot(3,2,1)
plot(t,u_L)
grid on
title('Left motor voltage')
xlabel('Time [s]')
ylabel('u_L [V]')

subplot(3,2,2)
plot(t,u_R)
grid on
title('Right motor voltage')
xlabel('Time [s]')
ylabel('u_R [V]')

subplot(3,2,3)
plot(t,v_L)
grid on
title('Left wheel velocity')
xlabel('Time [s]')
ylabel('v_L [m/s]')

subplot(3,2,4)
plot(t,v_R)
grid on
title('Right wheel velocity')
xlabel('Time [s]')
ylabel('v_R [m/s]')

subplot(3,2,5)
plot(t,volt)
grid on
title('Average voltage')
xlabel('Time [s]')
ylabel('volt [V]')

subplot(3,2,6)
plot(t,vel)
grid on
title('Average velocity')
xlabel('Time [s]')
ylabel('vel [m/s]')

%% =========================================================
% 7. IDENTIFY TRANSFER FUNCTIONS
% From voltage to wheel velocity
% ==========================================================

idd_L   = iddata(v_L,u_L,T_s);
idd_R   = iddata(v_R,u_R,T_s);
idd_avg = iddata(vel,volt,T_s);

% Tried 1 poles, 0 zero
G_L   = tfest(idd_L,1,0);
G_R   = tfest(idd_R,1,0);
G_avg = tfest(idd_avg,1,0);
[num_G_avg, den_G_avg] = tfdata(G_avg,'v');

%% =========================================================
% 8. COMPARE MODEL WITH MEASURED DATA
% ==========================================================

figure(101)

subplot(3,1,1)
compare(G_L,idd_L)
grid on
title('Compare - Left wheel')

subplot(3,1,2)
compare(G_R,idd_R)
grid on
title('Compare - Right wheel')

subplot(3,1,3)
compare(G_avg,idd_avg)
grid on
title('Compare - Average')

%% =========================================================
% 9. STEP RESPONSE: MODEL VS MEASURED DATA
% ==========================================================

[y_L,t_L]     = step(G_L,t(end));
[y_R,t_R]     = step(G_R,t(end));
[y_avg,t_avg] = step(G_avg,t(end));

figure(102)

subplot(3,1,1)
plot(t_L,y_L)
hold on
plot(t,v_L)
hold off
grid on
legend('Estimated output','Measured output')
title('Left wheel - Step response')
xlabel('Time [s]')
ylabel('Velocity [m/s]')

subplot(3,1,2)
plot(t_R,y_R)
hold on
plot(t,v_R)
hold off
grid on
legend('Estimated output','Measured output')
title('Right wheel - Step response')
xlabel('Time [s]')
ylabel('Velocity [m/s]')

subplot(3,1,3)
plot(t_avg,y_avg)
hold on
plot(t,vel)
hold off
grid on
legend('Estimated output','Measured output')
title('Average - Step response')
xlabel('Time [s]')
ylabel('Velocity [m/s]')

%% =========================================================
% 10. BODE PLOTS
% Week 6 - Compare wheels up / wheels down
% ==========================================================

figure(103)

subplot(1,2,1)
bode(G_avg)
grid on
title('Bode plot - wheels up')

% Data wheels down
% load and identify that transfer function as G_down
%
% Example:
% data2 = readtable('wheels_down.txt');
% ...
% G_down = tfest(idd_down,2,0);

% Temporary placeholder:
subplot(1,2,2)
%bode(G_down)
grid on
title('Bode plot - wheels down')

%% =========================================================
% 11. OPEN-LOOP WITH P-CONTROLLER
% G_ol = Kp * G
% ==========================================================

Kp = 15;     % must be <= 15

G_ol = Kp * G_avg;

figure(104)
margin(G_ol)
grid on
title('Open-loop with P-controller (Kp=15)')

%% =========================================================
% 12. LOW-PASS FILTER
% First-order filter:
% G_filt = 1 / (tau*s + 1)
% or
% G_filt = w_f / (s + w_f)
% ==========================================================
w1 = 10;
w2 = 100;

G_filt1 = tf(w1,[1 w1]);
G_filt2 = tf(w2,[1 w2]);

% Step response of filters
figure(105)

[y1,t1] = step(G_filt1);
[y2,t2] = step(G_filt2);

plot(t1,y1)
hold on
plot(t2,y2)
hold off
grid on

legend('w_b = 10','w_b = 100')
title('Low-pass filter step responses')
xlabel('Time [s]')
ylabel('Amplitude')



%% =========================================================
% 13. FILTERED OPEN-LOOP
% G_ol_filt = Kp * G * G_filt
% ==========================================================

G_ol_filt1 = Kp * G_avg * G_filt1;
G_ol_filt2 = Kp * G_avg * G_filt2;

figure(106)

subplot(2,1,1)
margin(G_ol_filt1)
grid on
title('Filtered open-loop (w_b = 10)')

subplot(2,1,2)
margin(G_ol_filt2)
grid on
title('Filtered open-loop (w_b = 100)')

% Step response of filtered open-loop
figure(107)

subplot(2,1,1)
step(G_ol_filt1)
grid on
title('Step response - filtered open-loop (w_b = 10)')

subplot(2,1,2)
step(G_ol_filt2)
grid on
title('Step response - filtered open-loop (w_b = 100)')

%% =========================================================
% 14. CLOSED-LOOP SYSTEM
% G_cl = (Kp*G) / (1 + Kp*G*G_filt)
% ==========================================================

G_cl = (Kp * G_avg) / (1 + Kp * G_avg * G_filt1);

figure(108)
step(G_cl)
grid on
title('Closed-loop step response')
xlabel('Time [s]')
ylabel('Amplitude')



%% =========================================================
% PI-Controller
% 15. Transfer function from voltage to position
% ==========================================================

G_pos = G_ol * tf(1, [1 0])

figure(109)
step(G_pos)
grid on
title('Voltage to position')
xlabel('Time [s]')
ylabel('Amplitude')


%% =========================================================
% 16. PI-LEAD CONTROL DESIGN FOR POSITION
% ==========================================================

figure(110)
bode(G_pos)
grid on
title('Bode plot - position system')

% Design parameters
N_i = 5;
alpha = 0.1;
gamma_M = 60;

% Frequency response
w = linspace(1e-2,120,1000);
[M,P,w_out] = bode(G_pos,w); 
M = mag2db(squeeze(M));
P = squeeze(P);

% Phase-balance equation
phi_i = rad2deg(-atan(1/N_i));
phi_G = -180 + gamma_M - phi_i;

% New crossover frequency
i_c = find(P <= phi_G,1,'first');
omega_c = w_out(i_c);

% PI part
tau_i = N_i/omega_c;
C_PI = tf([tau_i 1], [tau_i 0]);

% Lead part
tau_d = 1/(omega_c*sqrt(alpha));
C_D = tf([tau_d 1], [alpha*tau_d 1]);

% Open-loop without gain
G_pos_ol = minreal(C_PI * C_D * G_pos);

% Proportional gain
K_P = 1 / abs(squeeze(freqresp(G_pos_ol,omega_c)));

% Final open-loop
G_ol = minreal(K_P * C_PI * C_D * G_pos);

figure(111)
bode(G_ol)
grid on
title('Open-loop with PI-Lead controller')

%% =========================================================
% 17. CLOSED-LOOP - LEAD IN FORWARD BRANCH
% ==========================================================

L = K_P * C_PI * C_D * G_pos;
G_cl_a = L / (1 + L);

figure(112)
step(G_cl_a)

grid on
title('Closed-loop step response - Lead in forward branch')

%info_a = stepinfo(G_cl_a)
dcgain_a = dcgain(G_cl_a)

%% =========================================================
% 18. CLOSED-LOOP - LEAD IN FEEDBACK BRANCH
% ==========================================================

G_fw = K_P * C_PI * G_pos;
H_fb = C_D;
G_cl_b = G_fw / (1 + G_fw * H_fb);

figure(113)
step(G_cl_b)
grid on
title('Closed-loop step response - Lead in feedback branch')

%info_b = stepinfo(G_cl_b)
dcgain_b = dcgain(G_cl_b)

%% =========================================================
% 19. COMPARE STEP RESPONSES
% ==========================================================

figure(114)
step(G_cl_a,G_cl_b)
grid on
legend('Lead in forward branch','Lead in feedback branch')
title('Comparison of closed-loop step responses')
stepinfo;

%% =========================================================
% 20. Closed-loop step response from voltage to velocity
% ==========================================================

Kp = 4.68;
PI = Kp * tf(1, [1 0]);
OL = G_avg * PI;
CL = minreal(OL / (1+OL));

figure(115)
step(CL)
grid on
title('Closed-loop step response from voltage to velocity')






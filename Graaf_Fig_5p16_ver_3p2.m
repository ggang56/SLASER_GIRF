% Graaf_Fig_5p16_extra_ver_1p0.m
% Written by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: 17/09/2019, Last modified: 
%--------------------------------------------------------------------------
%% Clean slate
clear variables; close all; clc;

%% Add paths
addpath E:\mfiles_kmin\2DCAIPI\utils\export_fig

%% Set location
data_dir      = 'E:\mfiles_kmin\github\sLASER_GIRF-KMINNAM';
output_dir    = 'E:\mfiles_kmin\github\sLASER_GIRF-KMINNAM\output'; mkdir(output_dir);
data_filename = 'Graaf_Fig_5p16_ver_3p2';

% Robin_de_Graaf figure 5.16
dt = 6.4e-6; % RF dwell time  [sec]
T  = 4e-3;   % pulse duration [sec]

% Define discrete samples of time over [0 T)
nt = ceil(T/dt)+1;    % number of samples
t  = (0:nt-1).' * dt; % [sec]

% Set RF amplitude and frequency
b1max = 15; % [uT]

% define beta
beta = asech(0.01);    % [rad/msec]
mu   = 8.42;           % dimensionless parameter
vmax = mu*beta/(2*pi); % [kHz] mu*beta

% Figure 5.16 (A) and (B)
am_hs = b1max * sech(beta*(1-2*t/T)); % [uT]
fm_hs = vmax  * tanh(beta*(1-2*t/T)); % [kHz]

% Figure 5.16 (C)
% theta_t = cumtrapz(fm_hs); % degree of time
theta_t = cumtrapz(fm_hs) + 360/(2*pi)*T/2*mu*log(sech(beta*(1-2*t(1)/T))); % degree of time
% Figure 5.16 (D)
%--------------------------------------------------------------------------
% Original Catalina's script
% - fm_hs [Hz], dt [msec], beta [rad/msec], t [msec], T [msec]
% b_hs = am_hs.*exp(1j*(2*pi*cumtrapz(fm_hs)*dt*1e-3 + T/2*(2*pi*mu*1e-3)*log(sech(beta*(1-2*t(1)/T))))) .* exp(1j*2*pi*BW*offset*T); % [G]
%
% Modified 
% - define: fm_hs [kHz], dt [sec], beta [dimensionless], t [sec], T [sec]
% [1] remove the offset term and set the unit (uT) because of am_hs (uT)
% - 360 [degree] / [sec] =  2pi [rad] / [sec] = 1 [cycle] / [sec] = 1 [Hz]
% - cumtrapz(fm_hs) [degree][sec] -> 2*pi/360*cumtrapz(fm_hs) [radian][sec]
% b_hs = am_hs.*exp(1j*(2*pi/360*cumtrapz(fm_hs)+T/2*mu*log(sech(beta*(1-2*t(1)/T))))); % [uT]
% therefore, theta is cumtrapz(fm_hs) + 360/(2*pi)*T/2*mu*log(sech(beta*(1-2*t(1)/T)))); % degree of time
b_hs = am_hs .* exp(1j*(2*pi/360*theta_t)); % [uT]
b1x  = real(b_hs); b1y = imag(b_hs); % [uT]

%--------------------------------------------------------------------------
% Display
%--------------------------------------------------------------------------
LineWidth = 1.5;
black = [0 0 0]; gray  = [0.5 0.5 0.5];

% Figure 5.16 (A)
figure('color', 'w','units','normalized','outerposition',[0 0 1 1]); 
subplot(2,2,1); 
plot(t*1e3, am_hs, 'k', 'LineWidth', LineWidth);
xlabel('Time (ms)'); ylabel('RF amplitude (uT)');
set(gca, 'Box', 'off');

% Figure 5.16 (B)
subplot(2,2,2); 
plot(t*1e3, fm_hs,              'k',   'LineWidth', LineWidth); hold on;
plot(t*1e3, zeros(length(t),1), 'k--', 'LineWidth', LineWidth); hold off;
xlabel('Time (ms)'); ylabel('Frequency (kHz)'); set(gca, 'Box', 'off');

% Figure 5.16 (C)
subplot(2,2,3); 
plot(t*1e3, theta_t, 'k', 'LineWidth', LineWidth); hold on;
xlabel('Time (ms)'); ylabel('Phase (degree)'); set(gca, 'Box', 'off');

% Figure 5.16 (D)
subplot(2,2,4); 
plot(t*1e3, b1x,                'r',   'LineWidth', LineWidth); hold on;
plot(t*1e3, b1y,                'b',   'LineWidth', LineWidth); 
plot(t*1e3, zeros(length(t),1), 'k--', 'LineWidth', LineWidth); hold off;
xlabel('Time (ms)'); ylabel('RF amplitude (uT)'); set(gca, 'Box', 'off');
ylim([-15 15]); legend('B1x', 'B1y','Location','Southeast');

% save figure
export_fig(fullfile(output_dir, data_filename), '-m2', '-tif'); close;

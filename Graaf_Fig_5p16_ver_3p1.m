% Graaf_Fig_5p16_ver_3p1.m
% Written by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: 03/09/2019, Last modified: 17/09/2019
%--------------------------------------------------------------------------
%% Clean slate
clear variables; close all; clc;

%% Add paths
addpath E:\mfiles_kmin\2DCAIPI\utils\export_fig

%% Set location
data_dir      = 'E:\mfiles_kmin\github\sLASER_GIRF-KMINNAM';
output_dir    = 'E:\mfiles_kmin\github\sLASER_GIRF-KMINNAM\output'; mkdir(output_dir);
data_filename = 'Graaf_Fig_5p16_ver_3p1';

% Robin_de_Graaf figure 5.16
dt = 6.4e-6; % RF dwell time  [sec]
T  = 4e-3;   % pulse duration [sec]

% Define discrete samples of time over [0 T)
nt = ceil(T/dt)+1;    % number of samples
t  = (0:nt-1).' * dt; % [sec]

% Set RF amplitude and frequency
b1max = 15;           % [uT]
vmax  = 7.1;          % [kHz] mu*beta

% define beta
beta = asech(0.01);   % dimensionless parameter
mu   = 1;             % [kHz]

% Figure 5.16 (A) and (B)
am_hs = b1max * sech(beta*(1-2*t/T)); % [uT]
fm_hs = vmax  * tanh(beta*(1-2*t/T)); % [kHz]

% Figure 5.16 (C)
theta = cumtrapz(fm_hs) + 1/360*T/2*mu*log(sech(beta*(1-2*t(1)/T))); % degree

% Figure 5.16 (D)
b_hs = am_hs .* exp(1j*(2*pi/360*theta)); % [uT]
b1x  = real(b_hs); b1y = imag(b_hs);      % [uT]

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
plot(t*1e3, theta, 'k', 'LineWidth', LineWidth); hold on;
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






















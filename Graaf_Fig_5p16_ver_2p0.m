% Graaf_Fig_5p16_ver_2p0.m
% Written by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: 03/09/2019, Last modified: 17/09/2019
%--------------------------------------------------------------------------
clear variables; close all; clc;

data_dir      = 'F:\mfiles_kmin\github\sLASER_FOCI';
output_dir    = 'F:\mfiles_kmin\github\sLASER_FOCI\output'; mkdir(output_dir);
data_filename = 'Graaf_Fig_5p16_ver_2p0';

% Robin_de_Graaf figure 5.16
dt  = 6.4e-6;         % RF dwell time [sec]
tau = 4e-3;           % pulse duration, T [sec]

% Define discrete samples of time over [0 T)
nt = ceil(tau/dt)+1;  % number of samples
t  = (0:nt-1).' * dt; % [sec]

% Perform Bloch simulation
df = 1;         % sampling interval in freq. domain [Hz]
nf = 2e4;       % number of samples in freq. domain [-10, 10) kHz

% Set
b1max = 15;           % [uT]
vmax  = 7.1;          % [kHz] mu*beta

% define beta
beta = asech(0.01);

% Figure 5.16 (A) and (B)
AFP_am = b1max * sech(beta*(1-2*t/tau)); % [uT]
AFP_fm = vmax  * tanh(beta*(1-2*t/tau)); % [kHz]

% Figure 5.16 (C)
theta = cumtrapz(AFP_fm); % degree

% Figure 5.16 (D)
AFP_b1 = AFP_am .* exp(1i*2*pi/360*theta);  % [uT]
b1x = real(AFP_b1); b1y = imag(AFP_b1);     % [uT]

%--------------------------------------------------------------------------
% Display
%--------------------------------------------------------------------------
LineWidth = 1.5;
black = [0 0 0]; gray  = [0.5 0.5 0.5];

% Figure 5.16 (A)
figure('color', 'w','units','normalized','outerposition',[0 0 1 1]); 
subplot(2,2,1); plot(t*1e3, AFP_am, 'k', 'LineWidth', LineWidth);
xlabel('Time (ms)'); ylabel('RF amplitude (uT)');
set(gca, 'Box', 'off');

% Figure 5.16 (B)
subplot(2,2,2); plot(t*1e3, AFP_fm, 'k', 'LineWidth', LineWidth); hold on;
plot(t*1e3, zeros(length(t),1), 'k--', 'LineWidth', LineWidth); hold off;
xlabel('Time (ms)'); ylabel('Frequency (kHz)');
set(gca, 'Box', 'off');

% Figure 5.16 (C)
subplot(2,2,3); plot(t*1e3, theta, 'k', 'LineWidth', LineWidth);
xlabel('Time (ms)'); ylabel('Phase (degree)');
set(gca, 'Box', 'off');

% Figure 5.16 (D)
subplot(2,2,4); 
plot(t*1e3, b1x, 'r', 'LineWidth', LineWidth); hold on;
plot(t*1e3, b1y, 'b', 'LineWidth', LineWidth); 
plot(t*1e3, zeros(length(t),1), 'k--', 'LineWidth', LineWidth); hold off;
xlabel('Time (ms)'); ylabel('RF amplitude (uT)');
set(gca, 'Box', 'off'); legend('B1x', 'B1y','Location','Southeast');

save figure
export_fig(fullfile(output_dir, data_filename), '-m2', '-tif'); close;

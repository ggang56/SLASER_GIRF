% Graaf_Fig_5p16_ver_1p1.m
% Written by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: 03/09/2019, Last modified: 
%--------------------------------------------------------------------------
clear variables; close all; clc;

data_dir      = 'F:\mfiles_kmin\github\sLASER_FOCI';
output_dir    = 'F:\mfiles_kmin\github\sLASER_FOCI\output'; mkdir(output_dir);
data_filename = 'Graaf_Fig_5p16_ver_1p1';

% Robin_de_Graaf figure 5.16
dt  = 6.4e-6;         % RF dwell time [sec]
tau = 4e-3;           % pulse duration, T [sec]

% Define discrete samples of time over [0 T)
nt = ceil(tau/dt)+1;    % number of samples
t  = (0:nt-1).' * dt; % [sec]

% Perform Bloch simulation
df       = 1;   % sampling interval in frequency domain [Hz]
nf       = 2e4; % number of samples in frequency domain [-10, 10) kHz
Delta_Hz = (-floor(nf/2):ceil(nf/2)-1).' * df; % [Hz]

mx0 = zeros(nf,1, 'double');
my0 = zeros(nf,1, 'double');
mz0 = ones(nf,1, 'double');

% set 
b1max = 10;   % [uT]
vmax  = 7.1;  % [kHz]

% define beta
beta = asech(0.01);

fb = sech(beta*(1-2*t/tau));
fv = tanh(beta*(1-2*t/tau));

% Figure 5.16 (A) and (B)
AFP_am = b1max * fb;
AFP_fm = vmax * fv;

% Figure 5.16 (C)
theta = zeros(length(t),1);
for idx = 1:length(t)-1
   theta(idx+1,1) = theta(idx,1) + vmax * tanh(beta*(1-2*t(idx)/tau));
end

% Figure 5.16 (D)
b1 = AFP_am .* exp(1i*2*pi/360*theta);
b1x = real(b1); b1y = imag(b1);

% bloch simulation
[mx,my,mz] = bloch(b1*1e-2, zeros(nt,1), dt, 1e10, 1e10, Delta_Hz, 0, 0, mx0, my0, mz0);
mxy = mx + 1j*my;

% Display
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
plot(t*1e3, b1x, 'k', 'LineWidth', LineWidth); hold on;
plot(t*1e3, b1y, 'color', gray, 'LineWidth', LineWidth); hold off;
xlabel('Time (ms)'); ylabel('RF amplitude (uT)');
set(gca, 'Box', 'off'); legend('B1x', 'B1y','Location','Southeast');

% save figure
export_fig(fullfile(output_dir, data_filename), '-m2', '-tif'); close;

% Figure 5.16 extra
figure('color', 'w','units','normalized','outerposition',[0 0 1 1]); 
plot(Delta_Hz*1e-3, real(mxy), 'r', 'LineWidth', LineWidth); hold on;
plot(Delta_Hz*1e-3, imag(mxy), 'b', 'LineWidth', LineWidth);
plot(Delta_Hz*1e-3, abs(mxy), 'k', 'LineWidth', LineWidth); hold off;
xlabel('Frequency (kHz)'); ylabel('M_{xy}/M_0');
legend('M_x', 'M_y', 'M_{xy}', 'Location', 'Southeast');
set(gca, 'Box', 'off');

% save figure
export_fig(fullfile(output_dir, [data_filename, '_mxy']), '-m2', '-tif'); close;



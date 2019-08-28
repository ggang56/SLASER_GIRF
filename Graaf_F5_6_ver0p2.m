% Graaf_F5_6_ver0p1.m
% Written by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: 27/08/2019, Last modified: 
%--------------------------------------------------------------------------
clear variables; close all; clc;

data_dir      = 'F:\mfiles_kmin\sLASER_FOCI';
output_dir    = 'F:\mfiles_kmin\github\output';
data_filename = 'Graaf_figure_ch5p6';

% Robin_de_Graaf figure 5.6
dt  = 6.4e-6;   % RF dwell time [sec]
tau = 1e-3;     % pulse duration [sec]
n   = 3;        % zero crossing number

b1max = 0.18*1e2; % [uT], 1e4 G = 1T = 1e6 uT => 1G = 1e2 uT

% Define discrete samples of time over [0 T)
nt = ceil(tau/dt);    % number of samples
t  = (0:nt-1).' * dt; % [sec]

%% Perform Bloch simulation
df       = 1;   % sampling interval in frequency domain [Hz]
nf       = 2e4; % number of samples in frequency domain
Delta_Hz = (-floor(nf/2):ceil(nf/2)-1).' * df; % [Hz]
mx0 = zeros(nf,1); my0 = zeros(nf,1); mz0 = ones(nf,1);

fb_t      = sinc(2*n/tau*(t-(tau/2)));
b1_sinc_t = b1max * fb_t; % [uT]

% Calculate the excitation profile at on-resonance on SINC
% [mx,my,mz] = bloch(b1, gr, tp, t1, t2, df, dp, mode, mx, my, mz)
[mx_sinc,my_sinc,mz_sinc] = ...
    bloch(b1_sinc_t*1e-2, zeros(nt,1), dt, 1e2, 1e2, Delta_Hz, 0, 0, mx0, my0, mz0);

mxy_sinc = mx_sinc + 1j*my_sinc;
% scale_factor = max(abs(mxy_sinc));
% mxy_sinc = mxy_sinc / scale_factor;

%%
%--------------------------------------------------------------------------
% Display 5.6 (A, B, C)
%--------------------------------------------------------------------------
LineWidth = 1.5;
black = [0 0 0];
gray  = [0.5 0.5 0.5];

figure('color', 'w','units','normalized','outerposition',[0 0 1 1]); 
% Display 5.6 (A)
subplot(2,3,1), plot(t*1e3, b1_sinc_t, 'k', 'LineWidth', LineWidth); hold on;
plot(t*1e3, zeros(length(t),1), 'k--'); hold off; set(gca, 'Box', 'off');
xlabel('Time (ms)'); ylabel('RF amplitude (uT)');

% Display 5.6 (B)
subplot(2,3,2), hold on;
plot(Delta_Hz, mx_sinc, 'r', 'LineWidth', LineWidth);
plot(Delta_Hz, my_sinc, 'b', 'LineWidth', LineWidth);
plot(Delta_Hz, abs(mxy_sinc), 'k', 'LineWidth', LineWidth);
legend('M_x','M_y','M_{xy}','Location','Southeast'); hold off;
xlabel('Frequency (Hz)'); ylabel('M_{xy}/M_0');

% Display 5.6 (C)
subplot(2,3,3), plot(Delta_Hz, mz_sinc, 'k', 'LineWidth', LineWidth); hold on;
plot(Delta_Hz, zeros(length(Delta_Hz),1), 'k--'); hold off; set(gca, 'Box', 'off');
xlabel('Frequency (Hz)'); ylabel('M_z/M_0');

%%
% Robin_de_Graaf figure 5.6 (D)
beta = [2 7]; % b5 and b7
b1_gaussian_t = b1max * exp(-beta.*(2/tau*(t-tau/2)).^2);

% Calculate the excitation profile at on-resonance on Gaussian
% [mx,my,mz] = bloch(b1, gr, tp, t1, t2, df, dp, mode, mx, my, mz)
[mx_gaussian_beta_first,my_gaussian_beta_first,mz_gaussian_beta_first] = ...
    bloch(b1_gaussian_t(:,1)*1e-2, zeros(nt,1), dt, 1e10, 1e10, Delta_Hz, 0, 0, mx0, my0, mz0);
mxy_gaussian_beta_first = mx_gaussian_beta_first + 1j*my_gaussian_beta_first;

[mx_gaussian_beta_second,my_gaussian_beta_second,mz_gaussian_beta_second] = ...
    bloch(b1_gaussian_t(:,2)*1e-2, zeros(nt,1), dt, 1e10, 1e10, Delta_Hz, 0, 0, mx0, my0, mz0);
mxy_Gaussian_beta_second = mx_gaussian_beta_second + 1j*my_gaussian_beta_second;

%--------------------------------------------------------------------------
% Display 5.6 (D, E, F)
%--------------------------------------------------------------------------
% Display 5.6 (D)
subplot(2,3,4), plot(t, b1_gaussian_t(:,1), 'Color', black, 'LineWidth', LineWidth);
hold on; plot(t, b1_gaussian_t(:,2), 'Color', gray, 'LineWidth', LineWidth); hold off;
set(gca, 'Box', 'off'); xlabel('Time (ms)'); ylabel('RF amplitude');

% Display 5.6 (E)
subplot(2,3,5), hold on;
plot(Delta_Hz, mx_gaussian_beta_second, 'r', 'LineWidth', LineWidth);
plot(Delta_Hz, my_gaussian_beta_second, 'b', 'LineWidth', LineWidth);
plot(Delta_Hz, abs(mxy_Gaussian_beta_second), 'k', 'LineWidth', LineWidth);
legend('M_x','M_y','M_{xy}','Location','Southeast'); hold off;
xlabel('Frequency (kHz)'); ylabel('M_{xy}/M_0');

% Display 5.6 (F)
subplot(2,3,6), plot(Delta_Hz, mz_gaussian_beta_first, 'k', 'LineWidth', LineWidth); hold on;
plot(Delta_Hz, mz_gaussian_beta_second, 'Color', gray, 'LineWidth', LineWidth); 
plot(Delta_Hz, zeros(length(Delta_Hz),1), 'k--'); hold off; 
set(gca, 'Box', 'off'); xlabel('Frequency (kHz)'); ylabel('M_z/M_0');
legend('1% gauss','10% gauss','Location','Southeast');








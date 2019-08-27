% Graaf_F5_6_v001.m
% Written by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: 27/08/2019, Last modified: 
%--------------------------------------------------------------------------
clear variables; close all; clc;

data_dir      = 'F:\mfiles_kmin\sLASER_FOCI';
data_filename = 'Graaf_figure_ch5p6';

% Robin_de_Graaf figure 5.6
dt    = 6.4e-3; % [msec]
T     = 1;      % [msec]
n     = 3;      % zero crossing number
b1max = 0.15;   % [G] 1e4 G = 1T = 1e6 uT => 1G = 1e2 uT

% Define discrete samples of time over [0 T)
nt = floor(T/dt) + 1; % number of samples
T  = nt*dt;           % [msec] 
t  = (0:nt-1).' * dt; % [msec]

nvz = 512*2;
df = (-floor(nvz/2):ceil(nvz/2)-1).' / (nvz/2) * 10000;

%%
mx0 = zeros(nvz,1); 
my0 = zeros(nvz,1); 
mz0 = ones(nvz,1);

fb_t      = sinc(2*n/T*(t-T/2));
Sinc_b1_t = b1max * fb_t ;

slice_thickness = 4;  % slice thickness [cm]
BW = 1; % [Hz] * [kHz/1e3 Hz] => [kHz]
G  = BW / (4.257746778 * slice_thickness); % [kHz] / ([kHz/G] * [cm]) => [G/cm]
Sinc_g = G * ones(nt,1);

% Calculate the excitation profile at on-resonance on SINC
% [mx,my,mz] = bloch(b1, gr, tp, t1, t2, df, dp, mode, mx, my, mz)
[mx_Sinc_on_resonance,my_Sinc_on_resonance,mz_Sinc_on_resonance] = ...
    bloch(Sinc_b1_t, Sinc_g, dt*1e-3, 100, 100, df, 0, 0, mx0, my0, mz0);
mxy_Sinc_on_resonance = mx_Sinc_on_resonance + 1j*my_Sinc_on_resonance;

%--------------------------------------------------------------------------
% Display 5.6 (A, B, C)
%--------------------------------------------------------------------------
LineWidth = 1.5;
black = [0 0 0];
gray  = [0.5 0.5 0.5];

figure('color', 'w','units','normalized','outerposition',[0 0 1 1]); 
% Display 5.6 (A)
subplot(2,3,1), plot(t, Sinc_b1_t, 'k', 'LineWidth', LineWidth); hold on;
plot(t, zeros(length(t)), 'k--'); hold off; set(gca, 'Box', 'off');
xlabel('Time (ms)'); ylabel('RF amplitude (G)');

% Display 5.6 (B)
subplot(2,3,2), hold on;
plot(df*1e-3, mx_Sinc_on_resonance, 'r', 'LineWidth', LineWidth);
plot(df*1e-3, my_Sinc_on_resonance, 'b', 'LineWidth', LineWidth);
plot(df*1e-3, abs(mxy_Sinc_on_resonance), 'k', 'LineWidth', LineWidth);
legend('M_x','M_y','M_z','Location','Southeast'); hold off;
xlabel('Frequency (kHz)'); ylabel('M_{xy}/M_0');

% Display 5.6 (C)
subplot(2,3,3), plot(df*1e-3, mz_Sinc_on_resonance, 'k', 'LineWidth', LineWidth); hold on;
plot(df*1e-3, zeros(length(df)), 'k--'); hold off; set(gca, 'Box', 'off');
xlabel('Frequency (kHz)'); ylabel('M_z/M_0');

%%
% Robin_de_Graaf figure 5.6 (D)
beta = [5, 7]; % b5 and b7
Gaussian_b1_t = b1max * exp(-beta.*(2/T*(t-T/2)).^2);
Gaussian_g = G * ones(nt,1);

% Calculate the excitation profile at on-resonance on SINC
% [mx,my,mz] = bloch(b1, gr, tp, t1, t2, df, dp, mode, mx, my, mz)
[mx_Gaussian_on_resonance_b5,my_Gaussian_on_resonance_b5,mz_Gaussian_on_resonance_b5] = ...
    bloch(Gaussian_b1_t(:,1), Gaussian_g, dt*1e-3, 100, 100, df, 0, 0, mx0, my0, mz0);
mxy_Gaussian_on_resonance_b5 = mx_Gaussian_on_resonance_b5 + 1j*my_Gaussian_on_resonance_b5;

[mx_Gaussian_on_resonance_b7,my_Gaussian_on_resonance_b7,mz_Gaussian_on_resonance_b7] = ...
    bloch(Gaussian_b1_t(:,2), Gaussian_g, dt*1e-3, 100, 100, df, 0, 0, mx0, my0, mz0);
mxy_Gaussian_on_resonance_b7 = mx_Gaussian_on_resonance_b7 + 1j*my_Gaussian_on_resonance_b7;

%--------------------------------------------------------------------------
% Display 5.6 (D, E, F)
%--------------------------------------------------------------------------
% Display 5.6 (D)
subplot(2,3,4), plot(t, Gaussian_b1_t(:,1), 'Color', black, 'LineWidth', LineWidth);
hold on; plot(t, Gaussian_b1_t(:,2), 'Color', gray, 'LineWidth', LineWidth); hold off;
set(gca, 'Box', 'off');
xlabel('Time (ms)'); ylabel('RF amplitude');

% Display 5.6 (E)
subplot(2,3,5), hold on;
plot(df*1e-3, mx_Gaussian_on_resonance_b5, 'r', 'LineWidth', LineWidth);
plot(df*1e-3, my_Gaussian_on_resonance_b5, 'b', 'LineWidth', LineWidth);
plot(df*1e-3, abs(mxy_Gaussian_on_resonance_b5), 'k', 'LineWidth', LineWidth);
legend('M_x','M_y','M_z','Location','Southeast'); hold off;
xlabel('Frequency (kHz)'); ylabel('M_{xy}/M_0');

% Display 5.6 (F)
subplot(2,3,6), plot(df*1e-3, mz_Gaussian_on_resonance_b5, 'k', 'LineWidth', LineWidth); hold on;
plot(df*1e-3, mz_Gaussian_on_resonance_b7, 'Color', gray, 'LineWidth', LineWidth); 
plot(df*1e-3, zeros(length(df)), 'k--'); hold off; 
set(gca, 'Box', 'off'); xlabel('Frequency (kHz)'); ylabel('M_z/M_0');
legend('1% gauss','10% gauss','Location','Southeast');

%%
saveas(gcf, fullfile(data_dir, [data_filename,'.tif']));
close;







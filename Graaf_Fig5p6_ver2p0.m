% Graaf_F5_6_ver2p0.m
% Written by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: 27/08/2019, Last modified: 28/08/2019
%--------------------------------------------------------------------------
clear variables; close all; clc;

data_dir      = 'F:\mfiles_kmin\github\sLASER_FOCI';
output_dir    = 'F:\mfiles_kmin\github\sLASER_FOCI\output'; mkdir(output_dir);
data_filename = 'Graaf_figure_ch5p6';

% Robin_de_Graaf figure 5.6
dt  = 6.4e-6;   % RF dwell time [sec]
tau = 1e-3;     % pulse duration, T [sec]
n   = 3;        % zero crossing number

% Define discrete samples of time over [0 T)
nt = ceil(tau/dt);    % number of samples
t  = (0:nt-1).' * dt; % [sec]

alpha_excitation = 90;  % exciation flip angle [degrees]
alpha_inversion  = 180; % inversion flip angle [degrees]

fb = sinc(2*n/tau*(t-(tau/2)));
fb = fb / max(fb);

gamma = 4257.784679; % Gamma for current nucleus [Hz/G]
b1max_excitation = (alpha_excitation / 360) / (gamma * sum(fb) * dt) * 1e2; % [uT]
fprintf(sprintf('b1max_excitation = %f\n', b1max_excitation));

b1max_inversion = (alpha_inversion / 360) / (gamma * sum(fb) * dt) * 1e2; % [uT]
fprintf(sprintf('b1max_inversion = %f\n', b1max_inversion));

%% Perform Bloch simulation
df       = 1;   % sampling interval in frequency domain [Hz]
nf       = 2e4; % number of samples in frequency domain
Delta_Hz = (-floor(nf/2):ceil(nf/2)-1).' * df; % [Hz]
mx0 = zeros(nf,1); my0 = zeros(nf,1); mz0 = ones(nf,1);

b1_sinc_excitation = b1max_excitation * fb; % [uT]
b1_sinc_inversion = b1max_inversion * fb; % [uT]

% Calculate the excitation profile at on-resonance on SINC
% [mx,my,mz] = bloch(b1, gr, tp, t1, t2, df, dp, mode, mx, my, mz)
[mx_sinc_excitation,my_sinc_excitation,mz_sinc_excitation] = ...
    bloch(b1_sinc_excitation*1e-2, zeros(nt,1), dt, 1e2, 1e2, Delta_Hz, 0, 0, mx0, my0, mz0);

scale_factor = max(abs(mx_sinc_excitation));
mx_sinc_excitation = mx_sinc_excitation / scale_factor;

scale_factor = max(abs(my_sinc_excitation));
my_sinc_excitation = my_sinc_excitation / scale_factor;

mxy_sinc_excitation = mx_sinc_excitation + 1j*my_sinc_excitation;

[mx_sinc_inversion,my_sinc_inversion_,mz_sinc_inversion] = ...
    bloch(b1_sinc_inversion*1e-2, zeros(nt,1), dt, 1e2, 1e2, Delta_Hz, 0, 0, mx0, my0, mz0);

%% Display figures 5.6 (A, B, C)
LineWidth = 1.5;
black = [0 0 0];
gray  = [0.5 0.5 0.5];

figure('color', 'w','units','normalized','outerposition',[0 0 1 1]); 
% Display 5.6(A)
subplot(2,3,1), plot(t*1e3, b1_sinc_excitation, 'k', 'LineWidth', LineWidth); hold on;
plot(t*1e3, zeros(length(t),1), 'k--'); hold off; set(gca, 'Box', 'off');
xlabel('Time (ms)'); ylabel('RF amplitude (uT)');

% Display 5.6(B)
subplot(2,3,2), hold on;
plot(Delta_Hz, mx_sinc_excitation, 'r', 'LineWidth', LineWidth);
plot(Delta_Hz, my_sinc_excitation, 'b', 'LineWidth', LineWidth);
plot(Delta_Hz, abs(mxy_sinc_excitation), 'k', 'LineWidth', LineWidth);
legend('M_x','M_y','M_{xy}','Location','Southeast'); hold off;
axis([-1e4 1e4 -1.1 1.1]); xlabel('Frequency (Hz)'); ylabel('M_{xy}/M_0');

% Display 5.6(C)
subplot(2,3,3), plot(Delta_Hz, mz_sinc_inversion, 'k', 'LineWidth', LineWidth); hold on;
plot(Delta_Hz, zeros(length(Delta_Hz),1), 'k--'); hold off; set(gca, 'Box', 'off');
xlabel('Frequency (Hz)'); ylabel('M_z/M_0');

%%
% Robin_de_Graaf figure 5.6 (D)
beta1 = 4; % 1% gauss: the trancation of the asymptotic Gaussian funtion
beta2 = 2.5; % 10% gauss: the trancation of the asymptotic Gaussian funtion
gamma = 4257.784679; % Gamma for current nucleus [Hz/G]

gaussian_func = exp(-beta1.*(2/tau*(t-tau/2)).^2);

b1max_excitation = (alpha_excitation / 360) / (gamma * sum(gaussian_func) * dt) * 1e2; % [uT]
b1_gaussian_excitation_beta1 = b1max_excitation * exp(-beta1.*(2/tau*(t-tau/2)).^2);
b1_gaussian_excitation_beta2 = b1max_excitation * exp(-beta2.*(2/tau*(t-tau/2)).^2);

b1max_inversion = (alpha_inversion / 360) / (gamma * sum(gaussian_func) * dt) * 1e2; % [uT]
b1_gaussian_inversion_beta1 = b1max_inversion * exp(-beta1.*(2/tau*(t-tau/2)).^2);
b1_gaussian_inversion_beta2 = b1max_inversion * exp(-beta2.*(2/tau*(t-tau/2)).^2);

% Calculate the excitation profile at on-resonance on Gaussian
% [mx,my,mz] = bloch(b1, gr, tp, t1, t2, df, dp, mode, mx, my, mz)
% Excitation 
[mx_gaussian_excitation_beta1,my_gaussian_excitation_beta1,mz_gaussian_excitation_beta1] = ...
    bloch(b1_gaussian_excitation_beta1*1e-2, zeros(nt,1), dt, 1e10, 1e10, Delta_Hz, 0, 0, mx0, my0, mz0);
mxy_gaussian_excitation_beta1 = mx_gaussian_excitation_beta1 + 1j*my_gaussian_excitation_beta1;

[mx_gaussian_excitation_beta2,my_gaussian_excitation_beta2,mz_gaussian_excitation_beta2] = ...
    bloch(b1_gaussian_excitation_beta2*1e-2, zeros(nt,1), dt, 1e10, 1e10, Delta_Hz, 0, 0, mx0, my0, mz0);

scale_factor = max(abs(mx_gaussian_excitation_beta2));
mx_gaussian_excitation_beta2 = mx_gaussian_excitation_beta2 / scale_factor;

scale_factor = max(abs(my_gaussian_excitation_beta2));
my_gaussian_excitation_beta2 = my_gaussian_excitation_beta2 / scale_factor;

mxy_gaussian_excitation_beta2 = mx_gaussian_excitation_beta2 + 1j*my_gaussian_excitation_beta2;

% Inversion
[mx_gaussian_inversion_beta1,my_gaussian_inversion_beta1,mz_gaussian_inversion_beta1] = ...
    bloch(b1_gaussian_inversion_beta1*1e-2, zeros(nt,1), dt, 1e10, 1e10, Delta_Hz, 0, 0, mx0, my0, mz0);
mxy_gaussian_inversion_beta1 = mx_gaussian_inversion_beta1 + 1j*my_gaussian_inversion_beta1;

[mx_gaussian_inversion_beta2,my_gaussian_inversion_beta2,mz_gaussian_inversion_beta2] = ...
    bloch(b1_gaussian_inversion_beta2*1e-2, zeros(nt,1), dt, 1e10, 1e10, Delta_Hz, 0, 0, mx0, my0, mz0);

%--------------------------------------------------------------------------
% Display 5.6 (D, E, F)
%--------------------------------------------------------------------------
% Display 5.6 (D)
subplot(2,3,4), plot(t, b1_gaussian_excitation_beta1, 'Color', gray, 'LineWidth', LineWidth);
hold on; plot(t, b1_gaussian_excitation_beta2, 'Color', black, 'LineWidth', LineWidth); hold off;
set(gca, 'Box', 'off'); xlabel('Time (ms)'); ylabel('RF amplitude');

% Display 5.6 (E)
subplot(2,3,5), hold on;
plot(Delta_Hz, mx_gaussian_excitation_beta2, 'r', 'LineWidth', LineWidth);
plot(Delta_Hz, my_gaussian_excitation_beta2, 'b', 'LineWidth', LineWidth);
plot(Delta_Hz, abs(mxy_gaussian_excitation_beta2), 'k', 'LineWidth', LineWidth);
legend('M_x','M_y','M_{xy}','Location','Southeast'); hold off;
axis([-1e4 1e4 -1.1 1.1]); xlabel('Frequency (Hz)'); ylabel('M_{xy}/M_0');

% Display 5.6 (F)
subplot(2,3,6), plot(Delta_Hz, mz_gaussian_inversion_beta1, 'Color', gray, 'LineWidth', LineWidth); hold on;
plot(Delta_Hz, mz_gaussian_inversion_beta2, 'Color', black, 'LineWidth', LineWidth); 
plot(Delta_Hz, zeros(length(Delta_Hz),1), 'k--'); hold off; 
set(gca, 'Box', 'off'); xlabel('Frequency (Hz)'); ylabel('M_z/M_0');
legend('1% gauss','10% gauss','Location','Southeast');

saveas(gcf, fullfile(output_dir, [data_filename, '.tif']));
close;






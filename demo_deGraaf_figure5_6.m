% demo_deGraaf_figure5_6.m
% Written by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: 27/08/2019, Last modified: 
% Modified by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)

%% Clean slate
close all; clear all; clc;

%% Add paths
addpath D:\mfiles_nam\rf_pulse_design_code\Bloch_simulator;

%% Define parameters
dt = 6.4e-6;   % RF dwell time [sec]
T  = 1e-3;     % pulse duration [sec]
n  = 3;        % zero crossing number

% Define discrete samples of time over [0 T)
nt = ceil(T / dt);    % number of samples
t  = (0:nt-1).' * dt; % [sec]

gam  = 4257.784679 * 2 * pi * 1e-2; % [Hz/G] * [2*pi rad/cycle] * [G/1e2uT] => [rad/sec/uT]

%% Calculate an excitation sinc RF pulse
% [uT], 1e4 G = 1T = 1e6 uT => 1G = 1e2 uT
flipangle = 90 * pi / 180;  % flip angle [rad]
pulse  = sinc(2 * n / T * (t - (T / 2)));
b1max = flipangle / (gam * sum(pulse) * dt); % [uT]
b1_sinc_excitation = b1max * pulse; % [uT]

%% Calculate an inversion sinc RF pulse
% [uT], 1e4 G = 1T = 1e6 uT => 1G = 1e2 uT
flipangle = 180 * pi / 180;  % flip angle [rad]
pulse  = sinc(2 * n / T * (t - (T / 2)));
b1max = flipangle / (gam * sum(pulse) * dt); % [uT]
b1_sinc_inversion = b1max * pulse; % [uT]

%% Calculate excitation Gaussian RF pulses
cutoff1 = 1;  % [%]
cutoff2 = 10; % [%]
beta1 = -log(cutoff1 / 100); % truncated at an initial amplitude of 1%
beta2 = -log(cutoff2 / 100); % truncated at an initial amplitude of 10%

pulse1 = exp(-beta1 * (2 / T * (t - T / 2)).^2); % for t in [0,T)
pulse2 = exp(-beta2 * (2 / T * (t - T / 2)).^2); % for t in [0,T)

flipangle = 90 * pi / 180;  % flip angle [rad]
b1max1 = flipangle / (gam * sum(pulse1) * dt); % [uT]

b1_gaussian_excitation1 = b1max1 * pulse1; % [uT]
b1_gaussian_excitation2 = b1max1 * pulse2; % [uT]

%% Calculate inversion Gaussian RF pulses
% [uT], 1e4 G = 1T = 1e6 uT => 1G = 1e2 uT
flipangle = 180 * pi / 180;  % flip angle [rad]

b1max1 = flipangle / (gam * sum(pulse1) * dt); % [uT]
b1max2 = flipangle / (gam * sum(pulse2) * dt); % [uT]

b1_gaussian_inversion1 = b1max1 * pulse1; % [uT]
b1_gaussian_inversion2 = b1max2 * pulse2; % [uT]

%% Perform Bloch simulation
df       = 1;   % sampling interval in frequency domain [Hz]
nf       = 2e4; % number of samples in frequency domain [-10, 10) kHz
Delta_Hz = (-floor(nf/2):ceil(nf/2)-1).' * df; % [Hz]
mx0 = zeros(nf,1, 'double');
my0 = zeros(nf,1, 'double');
mz0 = ones(nf,1, 'double');

[mx,my,mz] = bloch(b1_sinc_excitation*1e-2, zeros(nt,1), dt, 1e10, 1e10, Delta_Hz, 0, 0, mx0, my0, mz0);
mxy_sinc_excitation = mx + 1j*my;

[mx,my,mz] = bloch(b1_sinc_inversion*1e-2, zeros(nt,1), dt, 1e10, 1e10, Delta_Hz, 0, 0, mx0, my0, mz0);
mz_sinc_inversion = mz;

[mx,my,mz] = bloch(b1_gaussian_excitation1*1e-2, zeros(nt,1), dt, 1e10, 1e10, Delta_Hz, 0, 0, mx0, my0, mz0);
mxy_gaussian_excitation1 = mx + 1j*my;

[mx,my,mz] = bloch(b1_gaussian_inversion1*1e-2, zeros(nt,1), dt, 1e10, 1e10, Delta_Hz, 0, 0, mx0, my0, mz0);
mz_gaussian_inversion1 = mz;

[mx,my,mz] = bloch(b1_gaussian_inversion2*1e-2, zeros(nt,1), dt, 1e10, 1e10, Delta_Hz, 0, 0, mx0, my0, mz0);
mz_gaussian_inversion2 = mz;

%% Display results
LineWidth = 1.5;
black = [0 0 0];
gray  = [0.5 0.5 0.5];
red  = [238 28  46 ] / 255;
blue = [0   83  159] / 255;

figure('color', 'w', 'Position', [1 81 1359 731]); 

% Display 5.6 (A)
subplot(2,3,1); hold on;
plot(t*1e3, b1_sinc_excitation, 'Color', black, 'LineWidth', LineWidth);
plot(t*1e3, zeros(length(t),1), 'k--');
xlabel('Time (ms)');
ylabel('RF amplitude (uT)');

% Display 5.6 (B)
subplot(2,3,2); hold on;
plot(Delta_Hz*1e-3, real(mxy_sinc_excitation), 'Color', red  , 'LineWidth', LineWidth);
plot(Delta_Hz*1e-3, imag(mxy_sinc_excitation), 'Color', blue , 'LineWidth', LineWidth);
plot(Delta_Hz*1e-3, abs(mxy_sinc_excitation) , 'Color', black, 'LineWidth', LineWidth);
legend('M_x', 'M_y', 'M_{xy}', 'Location', 'Southeast');
xlabel('Frequency (kHz)');
ylabel('M_{xy}/M_0');

% Display 5.6 (C)
subplot(2,3,3); hold on;
plot(Delta_Hz*1e-3, mz_sinc_inversion, 'Color', black, 'LineWidth', LineWidth);
plot(Delta_Hz*1e-3, zeros(length(Delta_Hz),1), 'k--');
xlabel('Frequency (kHz)');
ylabel('M_z/M_0');

% Display 5.6 (D)
subplot(2,3,4); hold on;
plot(t*1e3, b1_gaussian_excitation1, 'Color', gray , 'LineWidth', LineWidth);
plot(t*1e3, b1_gaussian_excitation2, 'Color', black, 'LineWidth', LineWidth);
xlabel('Time (ms)');
ylabel('RF amplitude (uT)');

% Display 5.6 (E)
subplot(2,3,5); hold on;
plot(Delta_Hz*1e-3, real(mxy_gaussian_excitation1), 'Color', red  , 'LineWidth', LineWidth);
plot(Delta_Hz*1e-3, imag(mxy_gaussian_excitation1), 'Color', blue , 'LineWidth', LineWidth);
plot(Delta_Hz*1e-3, abs(mxy_gaussian_excitation1) , 'Color', black, 'LineWidth', LineWidth);
xlabel('Frequency (kHz)');
ylabel('M_{xy}/M_0');

% Display 5.6 (F)
subplot(2,3,6); hold on;
plot(Delta_Hz*1e-3, mz_gaussian_inversion1, 'Color', gray , 'LineWidth', LineWidth);
plot(Delta_Hz*1e-3, mz_gaussian_inversion2, 'Color', black, 'LineWidth', LineWidth);
plot(Delta_Hz*1e-3, zeros(length(Delta_Hz),1), 'k--');
xlabel('Frequency (kHz)');
ylabel('M_z/M_0');
legend(sprintf('%d%% gauss', cutoff1), sprintf('%d%% gauss', cutoff2), 'Location', 'Southeast');

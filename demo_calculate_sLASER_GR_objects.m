% demo_calculate_sLASER_GR_objects.m
% Written by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/23/2019, Last modified: 08/23/2019

%% clean slate
close all; clear all; clc;

%% Define sLASER GR objects
demo_define_sLASER_GR_objects;

%% Define constants
MGG_GAMMA_1H = 42577.46778; % [Hz/mT]
DWELL        = 6.4e-3;      % gradient dwelltime [msec]

%% Calculate a hyperbolic secant pulse
%--------------------------------------------------------------------------
% Define parameters
%--------------------------------------------------------------------------
T     = 6.7;          % pulse duration [msec]
beta  = 1300;         % modulation angular frequency [rad/sec]
mu    = 4.9;          % dimensionless parameter
f     = 10;           % foci factor
bmax  = 20;           % [uT] 1e4G = 1T = 1e6uT => 1G = 1e2uT
gmax  = 40;           % maximum gradient strength in [mT/m] (usually 40 [mT/m] = 4 [G/cm])
smax  = 150;          % 20000 [G/cm/sec] * [1e1 mT/m / G/cm] * [sec/1e3 msec] => 200 [mT/m/msec]
slice_thickness = 20; % slice thickness [mm]
offset = 0;           % in units of slide widths

%--------------------------------------------------------------------------
% Define discrete samples of time over [-T/2 T/2)
%--------------------------------------------------------------------------
nt = floor(T / DWELL) + 1; % number of samples
T = nt * DWELL; % total time [msec]
t = (-floor(nt/2):ceil(nt/2)-1).' * DWELL; % [msec]

%--------------------------------------------------------------------------
% Calculate the bandwidth of a hyperbolic secant pulse
%--------------------------------------------------------------------------
bandwidth = 2 * mu * beta / (2*pi); % [dimensionless] * [rad/sec] / [rad] => [Hz]
gr_str = bandwidth / (MGG_GAMMA_1H * slice_thickness * 1e-3); % [Hz] / ([Hz/mT] * [mm] * [m/1e3 mm]) => [mT/m]

%--------------------------------------------------------------------------
% Calculate the frequency offset [Hz]
% Excitation k-space: image-domain (slice profile) <=> k-space (RF pulse)
% M(r - ro) <=> B1(t) * exp(-1j*k*ro), where k = -gamma * G * (T/2 - t) if t in [-T/2 T/2)
% => B1(t) * exp(-1j * -gamma * G * (T/2 - t) * zo), zo = n * dz
% =  B1(t) * exp( 1j *  gamma * G * dz * n * (T/2 - t)), BW = 1/(2*pi) * gamma * G * dz
% =  B1(t) * exp( 1j *  2 * pi * BW * n * (T/2 - t))
% =  B1(t) * exp( 1j * -2 * pi * BW * n * t) * exp( 1j *  2 * pi * BW * n * T/2)
%--------------------------------------------------------------------------
w_offset = -bandwidth * offset; % [Hz]

%--------------------------------------------------------------------------
% Calculate AM and FM shapes
%--------------------------------------------------------------------------
am_hs = bmax * sech(beta * t * 1e-3); % [uT]
fm_hs = -mu * beta / (2*pi) * tanh(beta * t * 1e-3) + w_offset; % [Hz]

%--------------------------------------------------------------------------
% Calculate the complex RF waveform
%--------------------------------------------------------------------------
b_hs = am_hs .* exp(1j * (2 * pi * cumtrapz(fm_hs) * DWELL * 1e-3 + mu * log(sech(beta * t(1) * 1e-3)) + 2 * pi * w_offset * t(1) * 1e-3)) * exp(1j * 2 * pi * bandwidth * offset * T/2 * 1e-3); % [uT]

%--------------------------------------------------------------------------
% Add slopes to RF and gradient waveforms
%--------------------------------------------------------------------------
slope_samples = ceil(gr_str / smax / DWELL); % [mT/m] / [mT/m/msec] / [msec]
leading = (0:slope_samples-1).'* DWELL * smax; % [msec] * [mT/m/msec] => [mT/m]
am_hs = [zeros(slope_samples,1); am_hs; zeros(slope_samples,1)]; % [uT]
fm_hs = [zeros(slope_samples,1); fm_hs; zeros(slope_samples,1)]; % [uT]
b_hs  = [zeros(slope_samples,1); b_hs;  zeros(slope_samples,1)]; % [uT]
g_hs  = [leading; gr_str * ones(nt,1); flip(leading)]; % [mT/m]

%%
lenc    = nt * DWELL;               % [msec]
slope   = slope_samples * DWELL;    % [msec]
dur     = lenc + 2 * slope;         % [msec]
time    = 100 * DWELL;              % [msec]
samples = nt + 2 * slope_samples;

s_echo = GR;
s_echo.dur         = dur;
s_echo.str         = gr_str;
s_echo.ori         = 0;       % M_ORI
s_echo.lenc        = lenc;
s_echo.slope       = slope;
s_echo.slope1      = slope;
s_echo.slope2      = slope;
s_echo.ref         = 0;       % reference point in the object, with respect to its begin [msec]
s_echo.time        = time;    % time of reference point of object within the sequence [msec]
s_echo.samples     = samples; % returns the number of samples that fit in the duration of the object
s_echo.gr_waveform = g_hs;

%% This code demonstrates how to add a gradient object to MPS gradient waveforms 
AQ_base_time = 35; % AQ`base:time: time of reference point of object within the sequence [msec]
% Note that the first RF starts before time "0". You need to take this into
% account when calculating "t".

nt = floor(AQ_base_time / DWELL);
t = (0:nt - 1).' * DWELL; % <- this part must be modified accordingly!

g1 = zeros(nt,1, 'double'); % M gradient waveform
g2 = zeros(nt,1, 'double'); % P gradient waveform
g3 = zeros(nt,1, 'double'); % S gradient waveform

am_shape = zeros(nt,1, 'double'); % amplitude-modulated waveform
fm_shape = zeros(nt,1, 'double'); % frequency-modulated waveform

% Attach s_echo to M gradient waveform
time_index = find(t == s_echo.time);
logical_index = false(nt,1);
logical_index((0:samples-1).' + time_index) = true;
g1(logical_index) = s_echo.gr_waveform;

% Attach am_hs to am_shape
am_shape(logical_index) = am_hs;

% Attach fm_hs to fm_shape
fm_shape(logical_index) = fm_hs;

%% Display gradients and RF waveforms
figure('Color', 'w');
subplot(5,1,1); plot(t, g1); grid on;
xlabel('Time (msec)');
ylabel('M (mT/m)');
ylim([0 gmax]);

subplot(5,1,4); plot(t, am_shape); grid on;
xlabel('Time (msec)');
ylabel('RF\_am (\muT)');
ylim([0 bmax]);

subplot(5,1,5); plot(t, fm_shape); grid on;
xlabel('Time (msec)');
ylabel('RF\_fm');






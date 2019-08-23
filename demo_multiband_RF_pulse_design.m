% demo_multiband_RF_pulse_design.m
% Written by Namgyun Lee
% Email: ggang56@gmail.com
% Created: 03/01/2017, Last modified: 03/01/2017

close all; clear all; clc;

%--------------------------------------------------------------------------
% Information about fremex05:
% amplitude range of am_fremex05: -163 <=> 32767 (from -32767 to 32767)
% amplitude range of fm_fremex05: -8192 <=> 5077 (from -8192 to 8192)
% number of samples: 500
%--------------------------------------------------------------------------
demo_read_RFVAR_archived_shapes;

%% Display the AM and FM of the high bandwidth "fremex05" single band RF pulse
am_shape = am_fremex05;
fm_shape = fm_fremex05;

% The contents of the float array <array-name> should be scaled between -1.0 and
% 1.0, so that at least one of these values occurs.
am_shape = am_shape / max(abs(am_shape));
fm_shape = fm_shape / max(abs(fm_shape));

dur = 6.5; % duration of the object (float [msec])
B1 = 20;   % maximum (circular) magnetic field strength of the pulse (float [uT])

% FM_scale: scale factor for the FM pulse shape (integer [Hz]). All elements in
% the FM pulse shape will be multiplied by this factor. The FM pulse shape
% itself is scaled between -1.0 and 1.0, with at least one of these values
% actually occurring.
FM_scale = 5e3; % (integer [Hz])
nr_samples = length(am_shape); % (= am_samples and fm_samples)

dt = dur / nr_samples; % [msec]
t = (0:nr_samples-1).' * dt; % [msec]

% Scale amplitude
am_shape = B1 * am_shape;
fm_shape = FM_scale * fm_shape;

figure;
subplot(2,2,1), plot(t, am_shape);
xlabel('time (ms)');
ylabel('B1 (\muT)');
axis([t(1) t(end) -B1 B1]);
title('single band');

subplot(2,2,3), plot(t, fm_shape*1e-3);
xlabel('time (ms)');
ylabel('frequency (kHz)');
axis([t(1) t(end) -FM_scale*1e-3 FM_scale*1e-3]);

%% Calculate the AM and FM of the high bandwidth "fremex05" dual band RF pulse
w0 = 5.5e3; % [Hz]

dt2 = dur / nr_samples * 2; % [msec]
t2 = (0:nr_samples-1).' * dt2; % [msec]

am_shape_dualband = am_shape .* cos(pi*w0*t2*1e-3);
fm_shape_dualband = fm_shape / 2 + w0;

subplot(2,2,2), plot(t2, am_shape_dualband);
xlabel('time (ms)');
ylabel('B1 (\muT)');
axis([t2(1) t2(end) -B1 B1]);
title('dual band');

subplot(2,2,4), plot(t2, fm_shape_dualband*1e-3);
xlabel('time (ms)');
ylabel('frequency (kHz)');
axis([t2(1) t2(end) 0 2*FM_scale*1e-3]);

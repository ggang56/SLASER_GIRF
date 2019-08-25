% demo_define_sLASER_GR_objects.m
% Written by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/22/2019, Last modified: 08/22/2019

close all; clear all; clc;


%% Define GR structure
GR = struct('dur'    , [], ... % duration of the object [msec]
            'str'    , [], ... % strength of the gradient [mT/m]
            'ori'    , [], ... % orientation of the gradients with respect to the physical gradient system.
                           ... % 'ori refers to one of the three orientation vectors, M, P, S. 
                           ... % Possible values: M_ORI, P_ORI, S_ORI
            'lenc'   , [], ... % length of the constant part. Only used if shape is TRAPEZOID [msec]
            'slope'  , [], ... % virtual attribute. Sets slope1 as well as slope2. Only used if shape is TRAPEZOID [msec] 
            'slope1' , [], ... % length of the leading slope. Only used if shape is TRAPEZOID [msec]
            'slope2' , [], ... % length of the trailing slope. Only used if shape is TRAPEZOID [msec]
            'ref'    , [], ... % reference point in the object, with respect to its begin [msec]
            'time'   , [], ... % time of reference point of object within the sequence [msec]
            'samples', []);    % returns the number of samples that fit in the duration of the object
                               % The value is calculated from the duration, samples = dur / dwell_dur + 2

%% Initialize all gradients
ex                  = GR;
r_ex                = GR;

%% Define gradient objects
% s_ex
s_ex.dur               = 7.3192;  % [msec]
s_ex.str               = -2.9699; % [mT/m]
s_ex.ori               = 0;       % M_ORI
s_ex.lenc              = 6.9192;  % [msec]
s_ex.slope             = 0.2000;  % [msec]
s_ex.slope1            = 0.2000;  % [msec]
s_ex.slope2            = 0.2000;  % [msec]
s_ex.ref               = 6.1563;  % [msec]
s_ex.time              = 0.0000;  % [msec]
s_ex.samples           = 1083;

laser_d_echo_spoil  = repmat(GR, [3 1]); % [0]: left  crusher on the 2nd FOCI pulse (M)
laser_r_echo_spoil  = repmat(GR, [3 1]); % [0]: right crusher on the 2nd FOCI pulse (M)
laser_d_echo        = GR;                %      left  crusher on the 2nd FOCI pulse (P)
laser_r_echo        = GR;                %      right crusher on the 2nd FOCI pulse (P) 

% d_echo_spoil1
d_echo_spoil1.dur      = 1.3000;   % [msec]
d_echo_spoil1.str      = -40.0000; % [mT/m]
d_echo_spoil1.ori      = 1;        % P_ORI
d_echo_spoil1.lenc     = 0.9000;   % [msec]
d_echo_spoil1.slope    = 0.2000;   % [msec]
d_echo_spoil1.slope1   = 0.2000;   % [msec]
d_echo_spoil1.slope2   = 0.2000;   % [msec]
d_echo_spoil1.ref      = 4.5832;   % [msec]
d_echo_spoil1.time     = 5.6165;   % [msec]
d_echo_spoil1.samples  = 142;

% s_echo
s_echo.dur             = 6.7328;   % [msec]
s_echo.str             = 0.0000;   % [mT/m]
s_echo.ori             = 1;        % P_ORI
s_echo.lenc            = 6.7328;   % [msec]
s_echo.slope           = 0.0000;   % [msec]
s_echo.slope1          = 0.0000;   % [msec]
s_echo.slope2          = 0.0000;   % [msec]
s_echo.ref             = 3.3536;   % [msec]
s_echo.time            = 5.6165;   % [msec]
s_echo.samples         = 1054;

% d_echo_spoil2
d_echo_spoil2.dur      = 1.3000;   % [msec]
d_echo_spoil2.str      = -40.0000; % [mT/m]
d_echo_spoil2.ori      = 2;        % S_ORI
d_echo_spoil2.lenc     = 0.9000;   % [msec]
d_echo_spoil2.slope    = 0.2000;   % [msec]
d_echo_spoil2.slope1   = 0.2000;   % [msec]
d_echo_spoil2.slope2   = 0.2000;   % [msec]
d_echo_spoil2.ref      = 4.5832;   % [msec]
d_echo_spoil2.time     = 14.1561;  % [msec]
d_echo_spoil2.samples  = 142;

% laser_s_echo
laser_s_echo.dur       = 6.7328;  % [msec]
laser_s_echo.str       = 0.0000;  % [mT/m]
laser_s_echo.ori       = 1;       % P_ORI
laser_s_echo.lenc      = 6.7328;  % [msec]
laser_s_echo.slope     = 0.0000;  % [msec]
laser_s_echo.slope1    = 0.0000;  % [msec]
laser_s_echo.slope2    = 0.0000;  % [msec]
laser_s_echo.ref       = 3.3536;  % [msec]
laser_s_echo.time      = 14.1561; % [msec]
laser_s_echo.samples   = 1054;

% laser_s_echo2
laser_s_echo2.dur      = 6.7328;  % [msec]
laser_s_echo2.str      = 0.0000;  % [mT/m]
laser_s_echo2.ori      = 2;       % S_ORI
laser_s_echo2.lenc     = 6.7328;  % [msec]
laser_s_echo2.slope    = 0.0000;  % [msec]
laser_s_echo2.slope1   = 0.0000;  % [msec]
laser_s_echo2.slope2   = 0.0000;  % [msec]
laser_s_echo2.ref      = 3.3536;  % [msec]
laser_s_echo2.time     = 30.6624; % [msec]
laser_s_echo2.samples  = 1054;

% r_echo_spoil2
r_echo_spoil2.dur      = 1.3000;  % [msec]
r_echo_spoil2.str      = 40.0000; % [mT/m]
r_echo_spoil2.ori      = 2;       % S_ORI
r_echo_spoil2.lenc     = 0.9000;  % [msec]
r_echo_spoil2.slope    = 0.2000;  % [msec]
r_echo_spoil2.slope1   = 0.2000;  % [msec]
r_echo_spoil2.slope2   = 0.2000;  % [msec]
r_echo_spoil2.ref      = -3.2832; % [msec]
r_echo_spoil2.time     = 21.2020; % [msec]
r_echo_spoil2.samples  = 142;

% d_echo2
d_echo2.dur            = 1.3000;  % [msec]
d_echo2.str            = 40.0000; % [mT/m]
d_echo2.ori            = 0;       % M_ORI
d_echo2.lenc           = 0.9000;  % [msec]
d_echo2.slope          = 0.2000;  % [msec]
d_echo2.slope1         = 0.2000;  % [msec]
d_echo2.slope2         = 0.2000;  % [msec]
d_echo2.ref            = 0.2000;  % [msec]
d_echo2.time           = 30.6624; % [msec]
d_echo2.samples        = 142;

% s_echo2
s_echo2.dur            = 6.7328;   % [msec]
s_echo2.str            = 0.0000;   % [mT/m]
s_echo2.ori            = 2;        % S_ORI
s_echo2.lenc           = 6.7328;   % [msec]
s_echo2.slope          = 0.0000;   % [msec]
s_echo2.slope1         = 0.0000;   % [msec]
s_echo2.slope2         = 0.0000;   % [msec]
s_echo2.ref            = 3.3536;   % [msec]
s_echo2.time           = 30.6624;  % [msec]
s_echo2.samples        = 1054;

% r_echo2
r_echo2.dur            = 1.3000;  % [msec]
r_echo2.str            = 40.0000; % [mT/m]
r_echo2.ori            = 0;       % M_ORI
r_echo2.lenc           = 0.9000;  % [msec]
r_echo2.slope          = 0.2000;  % [msec]
r_echo2.slope1         = 0.2000;  % [msec]
r_echo2.slope2         = 0.2000;  % [msec]
r_echo2.ref            = -3.2832; % [msec]
r_echo2.time           = 30.6624; % [msec]
r_echo2.samples        = 142;

% r_echo2_spoil1
r_echo2_spoil1.dur     = 1.3000;  % [msec]
r_echo2_spoil1.str     = 40.0000; % [mT/m]
r_echo2_spoil1.ori     = 1;       % P_ORI
r_echo2_spoil1.lenc    = 0.9000;  % [msec]
r_echo2_spoil1.slope   = 0.2000;  % [msec]
r_echo2_spoil1.slope1  = 0.2000;  % [msec]
r_echo2_spoil1.slope2  = 0.2000;  % [msec]
r_echo2_spoil1.ref     = -3.2832; % [msec]
r_echo2_spoil1.time    = 30.6624; % [msec]
r_echo2_spoil1.samples = 142;

%% Define constants
MGG_GAMMA_1H = 42577.46778; % [Hz/mT]
DWELL        = 6.4e-3;      % gradient dwelltime [msec]

%%
%--------------------------------------------------------------------------
%                   __________________
% GR`<name>:slope1 /                  \ GR`<name>:slope2
%                 /                    \ 
%                |------|---------------|
%               -3      0               5 ---> t (ms)
%                |------> GR`<name>:ref
%--------------------------------------------------------------------------
dur   = s_ex.dur;   % duration of the object [msec]
str   = s_ex.str;   % strength of the gradient [mT/m]
lenc  = s_ex.lenc;  % length of the constant part [msec]
slope = s_ex.slope; % virtual attribute. Sets slope1 as well as slope2 [msec]
ref   = s_ex.ref;   % reference point in the object, with respect to its begin [msec]
time  = s_ex.time;  % time of reference point of object within the sequence [msec]

%% Calculate read-only attributes
% duration of the second part of the object, between the reference point
% and the end [msec]
dur2 = dur - ref;

% zeroth moment of the first part of the gradient [mT/m * msec]
stat1 = str * (ref + ref - slope) / 2;

% zeroth moment of the second part of the gradient [mT/m * msec]
stat2 = str * (dur2 + dur2 - slope) / 2;

%--------------------------------------------------------------------------
% e.g.) slope = 3.7, DWELL = 1
%
%                  o---o---o---o---o---o
%               x  |   |   |   |   |   |  x
%           x   |  |   |   |   |   |   |  |   x
%       x   |   |  |   |   |   |   |   |  |   |   x
%   x---+---+---+--+-------------------+--+---+---+---x
%   0   1   2   3  3.7
%   1   2   3   4  1   2   3   4   5   6
%
% number of samples in slope = ceil(slope / DWELL) = 4
%
%
% e.g.) slope = 3, DWELL = 1
%
%               o-----------------o
%           x   |                     x
%       x   |   |                     |   x
%   x---+---+---+-----------------+---+---+---x
%   0   1   2   3
%   1   2   3
%
% number of samples in slope = ceil(slope / DWELL) = 3
%--------------------------------------------------------------------------
lenc_samples = lenc / DWELL

%(str * lenc) / DWELL





% stat1 = ref + lenc
% 
% 
% stat = str * (lenc + 2 * slope + lenc) / 2; % zeroth moment of the gradient [mT/m * msec]






% slope_samples = ceil(slope / DWELL);
% slope_updated = slope_samples * DWELL;
% 
% ref_samples = ceil((ref - slope + slope_updated) / DWELL);


%lenc_samples = ceil(lenc / DWELL)




% lenc = 
% 
% floor(s_ex.ref / DWELL) + 1
% 

% calculate the area and recalculate the gradient strength!


















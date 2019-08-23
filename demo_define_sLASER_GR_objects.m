% demo_define_sLASER_GR_objects.m
% Written by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/22/2019, Last modified: 08/22/2019

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

d_echo_spoil        = repmat(GR, [3 1]); % [0]: left  crusher on the 1st FOCI pulse (M)
r_echo_spoil        = repmat(GR, [3 1]); % [0]: right crusher on the 1st FOCI pulse (M)
d_echo              = GR;                %      left  crusher on the 1st FOCI pulse (P)
r_echo              = GR;                %      right crusher on the 1st FOCI pulse (P)

laser_d_echo_spoil  = repmat(GR, [3 1]); % [0]: left  crusher on the 2nd FOCI pulse (M)
laser_r_echo_spoil  = repmat(GR, [3 1]); % [0]: right crusher on the 2nd FOCI pulse (M)
laser_d_echo        = GR;                %      left  crusher on the 2nd FOCI pulse (P)
laser_r_echo        = GR;                %      right crusher on the 2nd FOCI pulse (P) 

d_echo2_spoil       = repmat(GR, [3 1]); % [1]: left  crusher on the 3rd FOCI pulse (P)
r_echo2_spoil       = repmat(GR, [3 1]); % [1]: right crusher on the 3rd FOCI pulse (P)
d_echo2             = GR;                %      left  crusher on the 3rd FOCI pulse (S)
r_echo2             = GR;                %      right crusher on the 3rd FOCI pulse (S)

laser_d_echo2_spoil = repmat(GR, [3 1]); % [1]: left  crusher on the 4th FOCI pulse (P)
laser_r_echo2_spoil = repmat(GR, [3 1]); % [1]: right crusher on the 4th FOCI pulse (P)
laser_d_echo2       = GR;                %      left  crusher on the 4th FOCI pulse (S)
laser_r_echo2       = GR;                %      right crusher on the 4th FOCI pulse (S)

%% Define gradient objects

% ex
ex.dur     = [];
ex.str     = [];
ex.ori     = 0; % M_ORI
ex.lenc    = [];
ex.slope   = [];
ex.slope1  = [];
ex.slope2  = [];
ex.ref     = [];
ex.time    = [];
ex.samples = [];

% r_ex
r_ex.dur     = [];
r_ex.str     = [];
r_ex.ori     = 0; % M_ORI
r_ex.lenc    = [];
r_ex.slope   = [];
r_ex.slope1  = [];
r_ex.slope2  = [];
r_ex.ref     = [];
r_ex.time    = [];
r_ex.samples = [];

% d_echo
d_echo.dur     = [];
d_echo.str     = [];
d_echo.ori     = 1; % P_ORI
d_echo.lenc    = [];
d_echo.slope   = [];
d_echo.slope1  = [];
d_echo.slope2  = [];
d_echo.ref     = [];
d_echo.time    = [];
d_echo.samples = [];

% r_echo
r_echo.dur     = [];
r_echo.str     = [];
r_echo.ori     = 1; % P_ORI
r_echo.lenc    = [];
r_echo.slope   = [];
r_echo.slope1  = [];
r_echo.slope2  = [];
r_echo.ref     = [];
r_echo.time    = [];
r_echo.samples = [];

% laser_d_echo
laser_d_echo.dur     = [];
laser_d_echo.str     = [];
laser_d_echo.ori     = 1; % P_ORI
laser_d_echo.lenc    = [];
laser_d_echo.slope   = [];
laser_d_echo.slope1  = [];
laser_d_echo.slope2  = [];
laser_d_echo.ref     = [];
laser_d_echo.time    = [];
laser_d_echo.samples = [];

% laser_r_echo
laser_r_echo.dur     = [];
laser_r_echo.str     = [];
laser_r_echo.ori     = 1; % P_ORI
laser_r_echo.lenc    = [];
laser_r_echo.slope   = [];
laser_r_echo.slope1  = [];
laser_r_echo.slope2  = [];
laser_r_echo.ref     = [];
laser_r_echo.time    = [];
laser_r_echo.samples = [];

% d_echo2
d_echo2.dur     = [];
d_echo2.str     = [];
d_echo2.ori     = 2; % S_ORI
d_echo2.lenc    = [];
d_echo2.slope   = [];
d_echo2.slope1  = [];
d_echo2.slope2  = [];
d_echo2.ref     = [];
d_echo2.time    = [];
d_echo2.samples = [];

% r_echo2
r_echo2.dur     = [];
r_echo2.str     = [];
r_echo2.ori     = 2; % S_ORI
r_echo2.lenc    = [];
r_echo2.slope   = [];
r_echo2.slope1  = [];
r_echo2.slope2  = [];
r_echo2.ref     = [];
r_echo2.time    = [];
r_echo2.samples = [];

% laser_d_echo2
laser_d_echo2.dur     = [];
laser_d_echo2.str     = [];
laser_d_echo2.ori     = 2; % S_ORI
laser_d_echo2.lenc    = [];
laser_d_echo2.slope   = [];
laser_d_echo2.slope1  = [];
laser_d_echo2.slope2  = [];
laser_d_echo2.ref     = [];
laser_d_echo2.time    = [];
laser_d_echo2.samples = [];

% laser_r_echo2
laser_r_echo2.dur     = [];
laser_r_echo2.str     = [];
laser_r_echo2.ori     = 2; % S_ORI
laser_r_echo2.lenc    = [];
laser_r_echo2.slope   = [];
laser_r_echo2.slope1  = [];
laser_r_echo2.slope2  = [];
laser_r_echo2.ref     = [];
laser_r_echo2.time    = [];
laser_r_echo2.samples = [];

for idx = 1:3
    ori = idx  - 1;

    % d_echo_spoil
    d_echo_spoil(idx).dur     = [];
    d_echo_spoil(idx).str     = [];
    d_echo_spoil(idx).ori     = ori;
    d_echo_spoil(idx).lenc    = [];
    d_echo_spoil(idx).slope   = [];
    d_echo_spoil(idx).slope1  = [];
    d_echo_spoil(idx).slope2  = [];
    d_echo_spoil(idx).ref     = [];
    d_echo_spoil(idx).time    = [];
    d_echo_spoil(idx).samples = [];

    % r_echo_spoil
    r_echo_spoil(idx).dur     = [];
    r_echo_spoil(idx).str     = [];
    r_echo_spoil(idx).ori     = ori;
    r_echo_spoil(idx).lenc    = [];
    r_echo_spoil(idx).slope   = [];
    r_echo_spoil(idx).slope1  = [];
    r_echo_spoil(idx).slope2  = [];
    r_echo_spoil(idx).ref     = [];
    r_echo_spoil(idx).time    = [];
    r_echo_spoil(idx).samples = [];
    
    % laser_d_echo_spoil
    laser_d_echo_spoil(idx).dur     = [];
    laser_d_echo_spoil(idx).str     = [];
    laser_d_echo_spoil(idx).ori     = ori;
    laser_d_echo_spoil(idx).lenc    = [];
    laser_d_echo_spoil(idx).slope   = [];
    laser_d_echo_spoil(idx).slope1  = [];
    laser_d_echo_spoil(idx).slope2  = [];
    laser_d_echo_spoil(idx).ref     = [];
    laser_d_echo_spoil(idx).time    = [];
    laser_d_echo_spoil(idx).samples = [];

    % laser_d_echo_spoil
    laser_r_echo_spoil(idx).dur     = [];
    laser_r_echo_spoil(idx).str     = [];
    laser_r_echo_spoil(idx).ori     = ori;
    laser_r_echo_spoil(idx).lenc    = [];
    laser_r_echo_spoil(idx).slope   = [];
    laser_r_echo_spoil(idx).slope1  = [];
    laser_r_echo_spoil(idx).slope2  = [];
    laser_r_echo_spoil(idx).ref     = [];
    laser_r_echo_spoil(idx).time    = [];
    laser_r_echo_spoil(idx).samples = [];

    % d_echo2_spoil
    d_echo2_spoil(idx).dur     = [];
    d_echo2_spoil(idx).str     = [];
    d_echo2_spoil(idx).ori     = ori;
    d_echo2_spoil(idx).lenc    = [];
    d_echo2_spoil(idx).slope   = [];
    d_echo2_spoil(idx).slope1  = [];
    d_echo2_spoil(idx).slope2  = [];
    d_echo2_spoil(idx).ref     = [];
    d_echo2_spoil(idx).time    = [];
    d_echo2_spoil(idx).samples = [];

    % r_echo2_spoil
    r_echo2_spoil(idx).dur     = [];
    r_echo2_spoil(idx).str     = [];
    r_echo2_spoil(idx).ori     = ori;
    r_echo2_spoil(idx).lenc    = [];
    r_echo2_spoil(idx).slope   = [];
    r_echo2_spoil(idx).slope1  = [];
    r_echo2_spoil(idx).slope2  = [];
    r_echo2_spoil(idx).ref     = [];
    r_echo2_spoil(idx).time    = [];
    r_echo2_spoil(idx).samples = [];

    % laser_d_echo2_spoil
    laser_d_echo2_spoil(idx).dur     = [];
    laser_d_echo2_spoil(idx).str     = [];
    laser_d_echo2_spoil(idx).ori     = ori;
    laser_d_echo2_spoil(idx).lenc    = [];
    laser_d_echo2_spoil(idx).slope   = [];
    laser_d_echo2_spoil(idx).slope1  = [];
    laser_d_echo2_spoil(idx).slope2  = [];
    laser_d_echo2_spoil(idx).ref     = [];
    laser_d_echo2_spoil(idx).time    = [];
    laser_d_echo2_spoil(idx).samples = [];

    % laser_r_echo2_spoil
    laser_r_echo2_spoil(idx).dur     = [];
    laser_r_echo2_spoil(idx).str     = [];
    laser_r_echo2_spoil(idx).ori     = ori;
    laser_r_echo2_spoil(idx).lenc    = [];
    laser_r_echo2_spoil(idx).slope   = [];
    laser_r_echo2_spoil(idx).slope1  = [];
    laser_r_echo2_spoil(idx).slope2  = [];
    laser_r_echo2_spoil(idx).ref     = [];
    laser_r_echo2_spoil(idx).time    = [];
    laser_r_echo2_spoil(idx).samples = [];
end

d_echo_spoil(1).str = 0;
d_echo_spoil(2).str = 0;
r_echo_spoil(1).str = 0;
r_echo_spoil(2).str = 0;

laser_d_echo_spoil(1).str = 0;
laser_d_echo_spoil(2).str = 0;
laser_r_echo_spoil(1).str = 0;
laser_r_echo_spoil(2).str = 0;

d_echo2_spoil(1).str = 0;
d_echo2_spoil(3).str = 0;
r_echo2_spoil(1).str = 0;
r_echo2_spoil(3).str = 0;

laser_d_echo2_spoil(1).str = 0;
laser_d_echo2_spoil(3).str = 0;
laser_r_echo2_spoil(1).str = 0;
laser_r_echo2_spoil(3).str = 0;

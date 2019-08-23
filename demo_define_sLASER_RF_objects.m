% demo_define_sLASER_RF_objects.m
% Written by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/22/2019, Last modified: 08/22/2019

%% Define RF structure
% You need to define it by yourself!
% Look at the maobjrfvar__mxg.c and figure out what parameters are needed
% to convert am_shape and fm_shape into amplitude in [uT] and pulse duration in [msec].
% Also look at how rfdtool converts "b1" to am_shape and fm_shape.
% You basically need to reverse this process.

%%
RF = struct('freq'    , [], ... % offset frequency with respect to the current F0 [Hz]
            'ref'     , [], ... % reference point in the object, with respect to its begin [msec]
            'time'    , [], ... % time of reference point of object within the sequence [msec]
            'dur'     , [], ... % duration of the object [msec]
                            ... % Note: dur is automatically rounded to multiples of 4 usec, the RF-dwell-time
            'phase'   , [], ... % offset phase with respect to resonance frequency F0 at the reference point [degrees]
            'AM_scale', [], ... % extra amplitude scale factor of the RF pulse
            'FM_scale', [], ... % scale factor for the FM pulse shape [Hz]. All elements in the FM pulse shape 
                            ... % will be multiplied by this factor. The FM pulse shape itself is scaled between -1.0 and 1.0, 
                            ... % with at least one of these values actually occurring.
            'angle'   , [], ... % pulse angle [degrees]
            'shape'   , [], ... % reference to an RFVAR object that defines the AM and/or FM waveforms
                         
        
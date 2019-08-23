% demo_read_RFVAR_archived_shapes.m
% Written by Namgyun Lee
% Email: ggang56@gmail.com
% Created: 03/01/2017, Last modified: 03/01/2017

close all; clear all; clc;

am_fremex05 = [
                 319,       358,       463,       543,       637, ...
                 721,       806,       877,       946,       995, ...
                1039,      1063,      1076,      1078,      1072, ...
                1059,      1049,      1043,      1052,      1070, ...
                1110,      1159,      1221,      1285,      1348, ...
                1405,      1451,      1483,      1500,      1502, ...
                1493,      1473,      1454,      1440,      1438, ...
                1460,      1502,      1571,      1653,      1749, ...
                1848,      1938,      2015,      2076,      2114, ...
                2132,      2126,      2108,      2079,      2049, ...
                2027,      2022,      2037,      2076,      2133, ...
                2205,      2279,      2345,      2403,      2436, ...
                2449,      2436,      2404,      2353,      2293, ...
                2232,      2186,      2153,      2152,      2173, ...
                2223,      2288,      2362,      2433,      2495, ...
                2541,      2564,      2568,      2549,      2514, ...
                2471,      2430,      2400,      2398,      2422, ...
                2488,      2586,      2709,      2847,      2986, ...
                3116,      3225,      3303,      3349,      3363, ...
                3348,      3315,      3277,      3250,      3244, ...
                3277,      3351,      3468,      3608,      3764, ...
                3909,      4037,      4126,      4174,      4173, ...
                4128,      4045,      3941,      3828,      3729, ...
                3665,      3644,      3675,      3750,      3859, ...
                3977,      4090,      4178,      4237,      4248, ...
                4222,      4154,      4060,      3951,      3844, ...
                3756,      3712,      3715,      3780,      3895, ...
                4055,      4241,      4431,      4612,      4768, ...
                4889,      4970,      5011,      5016,      5003, ...
                4978,      4974,      5002,      5088,      5232, ...
                5442,      5694,      5977,      6258,      6521, ...
                6740,      6910,      7018,      7071,      7079, ...
                7060,      7037,      7037,      7083,      7188, ...
                7360,      7583,      7843,      8111,      8369, ...
                8583,      8751,      8853,      8898,      8890, ...
                8851,      8803,      8773,      8789,      8874, ...
                9038,      9283,      9598,      9959,     10346, ...
               10729,     11093,     11414,     11695,     11925, ...
               12124,     12304,     12498,     12723,     13019, ...
               13391,     13867,     14419,     15058,     15730, ...
               16423,     17099,     17734,     18301,     18815, ...
               19250,     19654,     20042,     20427,     20896, ...
               21411,     22036,     22744,     23541,     24339, ...
               25171,     25938,     26654,     27233,     27762, ...
               28135,     28437,     28652,     28906,     29098, ...
               29420,     29738,     30181,     30612,     31214, ...
               31432,     32099,     32441,     32577,     32567, ...
               32548,     32445,     32281,     31935,     31636, ...
               31483,     31993,     31963,     31978,     32476, ...
               32602,     32628,     32592,     32479,     32313, ...
               31855,     31240,     30387,     29517,     28600, ...
               27792,     27102,     26635,     26340,     26240, ...
               26224,     26249,     26227,     26038,     25685, ...
               25098,     24283,     23259,     22102,     20844, ...
               19598,     18434,     17443,     16640,     16038, ...
               15605,     15250,     14897,     14469,     13904, ...
               13173,     12264,     11200,     10014,      8775, ...
                7547,      6422,      5489,      4822,      4438, ...
                4282,      4238,      4191,      4074,      3856, ...
                3546,      3192,      2882,      2729,      2813, ...
                3127,      3575,      4055,      4490,      4828, ...
                5043,      5128,      5098,      5009,      4944, ...
                5029,      5381,      6067,      7060,      8281, ...
                9637,     11023,     12370,     13584,     14629, ...
               15439,     16028,     16377,     16589,     16716, ...
               16947,     17389,     18225,     19453,     21110, ...
               23014,     25063,     27059,     28863,     30317, ...
               31414,     31924,     32025,     31723,     31251, ...
               30523,     29940,     29607,     29676,     30129, ...
               30886,     31827,     32673,     32761,     32767, ...
               32701,     32096,     30606,     28867,     26990, ...
               25312,     24145,     23741,     24220,     25379, ...
               27006,     28684,     30334,     31337,     32538, ...
               32645,     32549,     32269,     31999,     31680, ...
               31737,     31774,     32063,     32458,     32673, ...
               32724,     32695,     32456,     31462,     30304, ...
               28735,     27119,     25473,     23971,     22647, ...
               21526,     20497,     19476,     18276,     16817, ...
               15022,     12944,     10736,      8832,      8014, ...
                8991,     11368,     14374,     17458,     20310, ...
               22704,     24496,     25609,     25987,     25660, ...
               24691,     23205,     21397,     19501,     17818, ...
               16638,     16179,     16450,     17251,     18284, ...
               19265,     20001,     20345,     20254,     19717, ...
               18792,     17552,     16153,     14726,     13460, ...
               12503,     11972,     11848,     12034,     12360, ...
               12675,     12856,     12832,     12563,     12050, ...
               11316,     10401,      9365,      8277,      7230, ...
                6317,      5654,      5327,      5353,      5654, ...
                6107,      6599,      7047,      7397,      7620, ...
                7689,      7611,      7382,      7028,      6580, ...
                6098,      5642,      5298,      5137,      5193, ...
                5447,      5817,      6229,      6599,      6876, ...
                7016,      7004,      6830,      6511,      6066, ...
                5543,      4989,      4482,      4096,      3904, ...
                3923,      4112,      4391,      4683,      4927, ...
                5084,      5131,      5065,      4885,      4601, ...
                4234,      3798,      3323,      2825,      2329, ...
                1854,      1415,      1028,       693,       420, ...
                 203,        41,       -70,      -121,      -163
].';
             
fm_fremex05 = [
                -213,      -303,      -493,      -731,      -941, ...
               -1053,     -1132,     -1234,     -1404,     -1608, ...
               -1734,     -1750,     -1836,     -1999,     -2279, ...
               -2412,     -2591,     -2901,     -3198,     -3197, ...
               -3012,     -2819,     -2668,     -2489,     -2325, ...
               -2214,     -2153,     -2172,     -2260,     -2392, ...
               -2584,     -2821,     -3060,     -3276,     -3447, ...
               -3516,     -3459,     -3313,     -3111,     -2886, ...
               -2672,     -2469,     -2322,     -2257,     -2259, ...
               -2334,     -2475,     -2658,     -2866,     -3067, ...
               -3232,     -3335,     -3334,     -3257,     -3104, ...
               -2914,     -2714,     -2539,     -2396,     -2316, ...
               -2294,     -2327,     -2426,     -2566,     -2740, ...
               -2923,     -3077,     -3169,     -3186,     -3109, ...
               -2948,     -2739,     -2512,     -2312,     -2155, ...
               -2062,     -2029,     -2064,     -2155,     -2321, ...
               -2539,     -2784,     -2996,     -3136,     -3177, ...
               -3122,     -2980,     -2795,     -2594,     -2413, ...
               -2275,     -2203,     -2200,     -2258,     -2376, ...
               -2551,     -2770,     -3004,     -3220,     -3373, ...
               -3431,     -3375,     -3234,     -3028,     -2813, ...
               -2617,     -2466,     -2372,     -2346,     -2381, ...
               -2481,     -2633,     -2834,     -3055,     -3261, ...
               -3400,     -3437,     -3354,     -3170,     -2926, ...
               -2663,     -2427,     -2238,     -2114,     -2054, ...
               -2069,     -2144,     -2281,     -2461,     -2667, ...
               -2853,     -2986,     -3023,     -2955,     -2792, ...
               -2579,     -2359,     -2166,     -2027,     -1954, ...
               -1953,     -2014,     -2136,     -2307,     -2517, ...
               -2742,     -2960,     -3123,     -3209,     -3192, ...
               -3099,     -2936,     -2762,     -2592,     -2473, ...
               -2407,     -2412,     -2464,     -2560,     -2686, ...
               -2837,     -2996,     -3129,     -3212,     -3209, ...
               -3138,     -2992,     -2831,     -2652,     -2517, ...
               -2405,     -2360,     -2344,     -2388,     -2458, ...
               -2576,     -2700,     -2830,     -2915,     -2945, ...
               -2907,     -2814,     -2702,     -2564,     -2452, ...
               -2331,     -2274,     -2245,     -2294,     -2344, ...
               -2448,     -2549,     -2690,     -2791,     -2868, ...
               -2899,     -2884,     -2842,     -2767,     -2683, ...
               -2624,     -2560,     -2551,     -2576,     -2596, ...
               -2684,     -2759,     -2839,     -2914,     -2973, ...
               -3001,     -2992,     -2926,     -2883,     -2774, ...
               -2732,     -2650,     -2641,     -2582,     -2630, ...
               -2631,     -2706,     -2720,     -2823,     -2805, ...
               -2877,     -2831,     -2832,     -2737,     -2727, ...
               -2599,     -2604,     -2513,     -2538,     -2506, ...
               -2581,     -2589,     -2695,     -2719,     -2822, ...
               -2823,     -2856,     -2800,     -2787,     -2662, ...
               -2616,     -2525,     -2474,     -2433,     -2448, ...
               -2479,     -2551,     -2666,     -2805,     -2937, ...
               -3081,     -3157,     -3224,     -3156,     -3090, ...
               -2924,     -2810,     -2650,     -2538,     -2460, ...
               -2431,     -2465,     -2521,     -2664,     -2843, ...
               -3029,     -3225,     -3355,     -3385,     -3320, ...
               -3154,     -2922,     -2644,     -2425,     -2224, ...
               -2091,     -2046,     -2079,     -2235,     -2490, ...
               -2975,     -3636,     -4244,     -4609,     -4464, ...
               -3787,     -2811,     -1736,      -638,       510, ...
                1773,      3119,      4338,      5077,      5033, ...
                4232,      3061,      1952,      1074,       369, ...
                -263,      -927,     -1709,     -2602,     -3499, ...
               -4157,     -4272,     -3870,     -3198,     -2491, ...
               -2030,     -1679,     -1548,     -1498,     -1548, ...
               -1719,     -1920,     -2241,     -2636,     -3025, ...
               -3416,     -3654,     -3732,     -3548,     -3289, ...
               -2926,     -2656,     -2385,     -2265,     -2157, ...
               -2192,     -2221,     -2407,     -2595,     -2886, ...
               -3143,     -3425,     -3574,     -3629,     -3464, ...
               -3274,     -2952,     -2727,     -2485,     -2365, ...
               -2332,     -2377,     -2567,     -2841,     -3286, ...
               -3686,     -4029,     -4022,     -3654,     -2959, ...
               -2203,     -1468,      -822,      -224,       303, ...
                 810,      1279,      1723,      2106,      2414, ...
                2613,      2644,      2640,      2496,      2401, ...
                2258,      2209,      2180,      2225,      2300, ...
                2462,      2592,      2772,      2851,      2788, ...
                2632,      2288,      1848,      1259,       614, ...
                -238,     -1505,     -3658,     -6539,     -8192, ...
               -7207,     -4687,     -2795,     -1828,     -1383, ...
               -1159,     -1112,     -1097,     -1247,     -1415, ...
               -1719,     -2136,     -2671,     -3359,     -4142, ...
               -4704,     -4826,     -4517,     -3858,     -3216, ...
               -2674,     -2336,     -2086,     -2032,     -2017, ...
               -2183,     -2437,     -2798,     -3292,     -3797, ...
               -4184,     -4259,     -3989,     -3478,     -2885, ...
               -2349,     -1942,     -1649,     -1447,     -1367, ...
               -1377,     -1468,     -1755,     -2204,     -2878, ...
               -3704,     -4412,     -4628,     -4265,     -3531, ...
               -2695,     -1972,     -1502,     -1230,     -1092, ...
               -1068,     -1173,     -1390,     -1713,     -2171, ...
               -2786,     -3475,     -4054,     -4347,     -4236, ...
               -3798,     -3209,     -2660,     -2231,     -1943, ...
               -1774,     -1738,     -1804,     -1966,     -2254, ...
               -2762,     -3456,     -4164,     -4749,     -5004, ...
               -4739,     -4086,     -3325,     -2607,     -2024, ...
               -1607,     -1354,     -1180,      -990,      -811, ...
                -841,      -880,      -758,      -890,     -1032, ...
                -943,      -979,      -735,      -388,        98, ...
                 513,       788,       971,      1046,      1030 ...
].';

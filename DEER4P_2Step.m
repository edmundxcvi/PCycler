% Sample input for pcycler
% 4PDEER experiment with conventional 2-step phase cycle
clear all

% Experiment variables - REQUIRED
% Number of pulses
Exp.nPulses = 4;
% Desired coherence transfer pathway
% n.b. Starts at equilibrium and ends at detection
Exp.p_des = [0 -1 +1 +1 -1];
% Names for interpulse delays
syms tau1 tau2 t_DEER
Exp.d_names = [tau1, tau1 + t_DEER, tau2 - t_DEER];
Exp.var = t_DEER;

% Experiment variables - OPTIONAL
% Phase cycle
Exp.p_cycle = [0      0      0      0      0      ;
               pi     0      0      0      0      ];
% Maximum coherence order
Exp.max_order = 1;
% Initial coherence order
Exp.p0 = 0;

% Simulation options - OPTIONAL
% Ingore FIDs
Opt.FID = 0;
% Include echos which do not cross detection
Opt.onlycrossing = 1;

% Coherence filtering example
% n.b. this is unnecessary and included for demonstration purposes only.
% Detection coherence order is -1 automatically.
% n.b. Starts at equilibrium and ends at detection
Opt.filter = [-1i, -1i, -1i, -1];


% Run calculation
[ctp, cycled_weights] = pcycler(Exp, Opt);

       
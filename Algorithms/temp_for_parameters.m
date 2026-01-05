k = 4;
% BH IC: -----------------------------------------------
S0_values = [100, 200, 400, 800];
I0_values = [10, 20, 40, 80];
S0 = S0_values(k); I0 = I0_values(k);
SI0 = [S0; I0]; 

% BH parameters: -----------------------------------------------
BH_parms = struct;
BH_parms.mu = 10 ^ (-5);           % natural death rate
BH_parms.lambda = 10 ^ (-5);       % S introduction rate
% Coefficient linking the WH viral load to the individual disease
% transmission rate:
l_values = [10^(-9.8), 10^(-9.8), 10^(-9.8), 10^(-9.8)]; 
BH_parms.l = l_values(k);

% WH IC: -----------------------------------------------
T0_infect = 10 ^ (8);          % healthy target cells
Tstar0_infect = 10 ^ (0);      % infectious cells
V0_infect = 10 ^ (6);          % viral particles
TIV0_infect = [T0_infect; Tstar0_infect; V0_infect]; % this are initial
                               % conditions for the INFECTIOUS agents

% WH parameters: -----------------------------------------------
WH_parms = struct;
WH_parms.p = 10 ^ (1);           % viral production rate (p per day)
WH_parms.c = 10 ^ (0);           % viral clearance rate (c per day)
WH_parms.k = 10 ^ (-7);          % infection coefficient of target cells 
                                 % (k per virus per day)
WH_parms.mu_c = 10 ^ (-1);       % mortality rate of cells (mu_c per day)
WH_parms.delta_c = 10 ^ (1);     % extra mortality rate of infectious cells 
                                 % (delta_c per day)
WH_parms.R_0_WH = 1.1;           % within-host basic reproduction number
WH_parms.Lambda_c = WH_parms.R_0_WH * ...
    (WH_parms.mu_c * (WH_parms.mu_c + WH_parms.delta_c) * WH_parms.c)...
    / (WH_parms.k * WH_parms.p);

% Simulation setup: -----------------------------------------------
t0 = 0; 
t_end_values = [20, 20, 20, 20];
t_endBH = t_end_values(k);
t_endWH = t_endBH; 

% Resolution to extract the poplation-level SI dynamics over 800
% repetitions:
rep = 800; 
delta_t_uniform = 0.1;
t_uniform = (t0 : delta_t_uniform : t_endBH);

% The resolution to harvest the within-host information: ------------------
pt_N = 130;
delta_t_values = [0.0000937 * 1.084.^(0:pt_N);
    0.0000937 * 1.084.^(0:pt_N);
    0.0000937 * 1.084.^(0:pt_N);
    0.0000937 * 1.084.^(0:pt_N)];
% Resolution to extract WH info: ------------------------------------------
delta_t = delta_t_values(k, :);     

[t_WH, V_infect, Psi_table, const_Psi, t_table] = WH_exact(TIV0_infect, WH_parms, t_endWH, 0.0000937, BH_parms);

BH_parms.l*max(V_infect)*((880)^2)/4 + BH_parms.lambda + BH_parms.mu*800 + BH_parms.mu*80
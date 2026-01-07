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

[t_WH, V_infect, C_infect, Cstar_infect, Psi_table, const_Psi, t_table] = ...
    WH_exact(TIV0_infect, WH_parms, t_endWH, 0.0000937, BH_parms);

BH_parms.l*max(V_infect)*((880)^2)/4 + BH_parms.lambda + BH_parms.mu*800 + BH_parms.mu*80

% plot the within-host dynamics
figure;
set(gcf, 'Position', get(0, 'Screensize'));
    ax = gca;
    ax.FontSize = 12;      
    ax.LabelFontSizeMultiplier = 1;
    ax.TitleFontSizeMultiplier = 1.2;
set(groot,'defaultAxesFontName','Verdana');
set(groot,'defaultAxesFontSize',20);
set(0, 'DefaultLineLineWidth', 4);
set(groot,'defaultLineMarkerSize', 6);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(t_WH, C_infect, 'DisplayName','$C$ (cell num.)'); hold on;
plot(t_WH, Cstar_infect, 'DisplayName','$C^*$ (infected cell num.)');
plot(t_WH, V_infect, 'DisplayName','$V$ (viral part. num.)');
xlabel('infection age $\delta t$ (days)');
ylabel('within-host quantities');
xlim([0 10]);
%title('Viral Load Dynamics Over Time');
legend('Location','best');
    lgd.FontSize = 14;
grid on;
hold off;


function dydt = WH_infect_odesystem_Log(t, y, WH_parms)
% The ODE equations describing the dynamics of WH subsystem among I populations. 
% Note that:
% y(1) := T
% y(2) := Tstar
% y(3) := V
dydt = [
    WH_parms.Lambda_c*exp(-y(1)) - WH_parms.k*exp(y(3)) - WH_parms.mu_c;
    WH_parms.k*exp(y(1)+y(3)-y(2)) - WH_parms.mu_c - WH_parms.delta_c;
    -WH_parms.c + WH_parms.p*exp(y(2)-y(3))
];
end

function [t_WH, V_infect, C_infect, Cstar_infect, Psi_table, const_Psi, t_table] = ...
    WH_exact(TIV0_infect, WH_parms, t_endWH, delta_t, BH_parms)
% This function solves the TIV model for the infectious population and
% outputs V_infect depending on the age of infection from 0 to 't_endWH'.

Initial_Conds_infect = log(TIV0_infect); % log trandsform the ICs to solve ODE
% time points based on uniform t discretisation for extracting WH info:
t_WH = (delta_t/2 : delta_t : ceil(2*t_endWH/delta_t)*delta_t/2);
t_WH = [0, t_WH]; % include ICs to solve the WH ODE.

opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[~, y1] = ode15s(@(t,y) WH_infect_odesystem_Log(t, y, WH_parms), t_WH, ...
    Initial_Conds_infect, opts);
t_WH = t_WH(2 : end);
V_infect_Log = y1(: , 3); 
V_infect = exp(V_infect_Log);
V_infect = V_infect(2 : end); % WH viral load at t_WH

C_infect_Log = y1(: , 1); 
C_infect = exp(C_infect_Log);
C_infect = C_infect(2 : end); % WH healthy cells at t_WH


Cstar_infect_Log = y1(: , 2); 
Cstar_infect = exp(Cstar_infect_Log);
Cstar_infect = Cstar_infect(2 : end); % WH infected cells at t_WH



% Record the 'Psi' table based on uniform t discretisation:
Psi_table = cumtrapz(t_WH, V_infect) * BH_parms.l; % approximate cumulative integral

% Record 'Psi' with constant probability discretisation but variable t
% values:
Psi_table_min = min(Psi_table);
Psi_table_max = max(Psi_table);
const_Psi = linspace(Psi_table_min, Psi_table_max, length(Psi_table));  
% obtatain variable time points at which WH info is extracted:
t_table = interp1(Psi_table, t_WH, const_Psi, 'linear', 'extrap');

end
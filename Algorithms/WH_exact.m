function [t_WH, V_infect, Psi_table, const_Psi, t_table] = ...
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
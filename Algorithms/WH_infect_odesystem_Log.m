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
function [T, Dynamics] = SI_FRM(initialState, tMax, dT, Lambda, beta, mu)

% Initialize time vector and dynamics matrix
time = 0:dT:tMax;
Dynamics = zeros(length(time), length(initialState));
T = time;

% Set the initial state in the dynamics matrix
T(1) = 0;
Dynamics(1, :) = initialState;
currentTime = 0;
currentState = initialState;

% Computation of reaction propensities:
a = zeros(1,4);
a(1) = Lambda;
a(2) = mu*currentState(1);
a(3) = mu*currentState(2);
a(4) = beta*currentState(1)*currentState(2);
a0 = sum(a);

% Iterate over time steps to compute dynamics
for i = 2:length(T)
    currentTime = currentTime + dT;
    
    if (rand <= a0*dT)

        % Update the current state based on the selected reaction
        reactionIndex = find(rand <= cumsum(a)/a0, 1);
        switch reactionIndex
            case 1
                currentState(1) = currentState(1) + 1; % Reaction 1
            case 2
                currentState(1) = currentState(1) - 1; % Reaction 2
            case 3
                currentState(2) = currentState(2) - 1; % Reaction 3
            case 4
                currentState(1) = currentState(1) - 1; 
                currentState(2) = currentState(2) + 1; % Reaction 4
        end
        
        % Update the total propensity for the next iteration
        a(2) = mu * currentState(1);
        a(3) = mu * currentState(2);
        a(4) = beta * currentState(1) * currentState(2);
        a0 = sum(a);
        Dynamics(i, :) = currentState; % Store the current state in the dynamics matrix
    else
    Dynamics(i,:) = Dynamics(i-1,:);
    end
end
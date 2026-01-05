for i=1:50

% initial conditions and initialisations
i0=20;    % initial condiction for I
s0=80; % initial condiction for S
tMax=50;   % evaluation time
dT= 0.1;

% parameters 
Lambda = 10^-1; % birth rate
beta=1.5*10^-3; % infectious rate
mu = 10^-5; % death rate

initialState=[s0 i0];

% FRM algorithm
[T, SI_FRM_Dynamics] = SI_FRM(initialState, tMax, dT, Lambda, beta, mu);

plot(T,SI_FRM_Dynamics(:,1),'b');
hold on;
grid on;
plot(T,SI_FRM_Dynamics(:,2),'r');
title('SI stochastic model');
end

% Numerical Integration
[T,Y] = ode45(@(t,Y) SIequations(t,Y,Lambda, beta, mu),Tspam,S0I0);

S=Y(:,1); % Solution S
I=Y(:,2); % Solution I

plot(T,S,'k', 'LineWidth',2);
hold on;
grid on;
plot(T,I,'k', 'LineWidth',2);
title(['SI model with parameters: \beta= ',num2str(beta),' N=',num2str(N)])
xlabel('Time')
ylabel('Number of Individuals')
legend('S','I','Location','best')
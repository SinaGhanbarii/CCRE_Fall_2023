%
%--CCRE (Prof. Matteo Maestri) - 2023/2024 - October 23  
%
% reaction kinetics:
% A -> P -> C
%
% variables: 
% j = 1 is A || j = 2 is P || j = 3 is C || j = 4 is T

clear all, clc

global E1 R k1_0 E2 k2_0 deltaH1 deltaH2 U A_v rho cp t_s

scale = 'industrial'; % industrial || lab

%% Data
rho = 1000; % kg/m3
cp = 4*1e3; %J/kg/K
U = 500; %W/m2/K
k1_0=0.5; % 1/s
k2_0=1e11; % 1/s
E1=20e3; %J/mol
E2=100e3; %J/mol
deltaH1 = -300e3; %J/mol
deltaH2 = -250e3; 
R = 8.314; % J/mol/K
t_s = 3600; % initial value: 3600 s

switch scale
    case 'lab'
    % lab scale values:
    A_v = 9.5; % 1/m
    V = 0.063; % m3
    case 'industrial'
    % industrial scale values:
    A_v = 2.6; % 1/m
    V = 6.3; % m3
end

% initial conditions:
T0 = 295; % K
cA_0=1000; % mol/m3
cP_0 = 0;
cC_0 = 0;

tspan=[0 6000]; % seconds
[t,y] = ode23s(@batch_reactor, tspan, [cA_0 cP_0 cC_0 T0]);

% post-processing:
cA = y(:,1);
cP = y(:,2);
cC = y(:,3);
T = y(:,4);

T_h = zeros(length(t),1);

for i = 1:length(t)
    if t(i)<t_s
        T_h(i) = 345; %K
    else 
        T_h(i) = 295; %K
    end
end

Conversion = 1-cA(end)/cA_0;
Selectivity = cP(end)/(cA_0-cA(end));
Yield = Conversion*Selectivity;

% Temperature profiles
figure(1)
plot(t,T)
hold on
xlabel('Time [s]','FontSize',14)
ylabel('Reactor temperature [K]','FontSize',14)
legend(scale)

tt = [t_s, t_s];
Treact = [min(T), max(T)];

plot(tt, Treact, '--')
% concentration profiles
switch scale
    case 'lab'
    figure(2)
    plot(t,y(:,1:3))
    hold on
    xlabel('Time [s]','FontSize',14)
    ylabel('Concentration [mol/m3]','FontSize',14)
    legend('A','P','C')
    tt = [t_s, t_s];
    creact = [0,max(cA)];
    plot(tt, creact, '--')

    case 'industrial'
    figure(3)
    plot(t,y(:,1:3))
    hold on
    xlabel('Time [s]','FontSize',14)
    ylabel('Concentration [mol/m3]','FontSize',14)
    legend('A','P','C')
    tt = [t_s, t_s];
    creact = [0,max(cA)];
    plot(tt, creact, '--')

end

figure(4)
plot(t,cP)
hold on
xlabel('Time [s]','FontSize',14)
ylabel('Concentration of P [mol/m3]','FontSize',14)
ylim([0 600])
legend(scale)

fprintf("Conversion: %f %% \n",Conversion*100)
fprintf("Selectivity: %f %% \n",Selectivity*100)
fprintf("Yield: %f %% \n",Yield*100)


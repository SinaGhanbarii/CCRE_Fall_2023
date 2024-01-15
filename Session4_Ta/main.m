close, clear, clc
global V c0 c_0 Q rho cp delH_rxn T_feed nu UA k0 Tcool

% Initial data
V = 10; % m^3
c0 = 5; % kmol/m^3
Q = 0.01; % m^3/s
k0 = 10^13; % 1/s
rho = 850; % kg/m^3
cp = 2200; % kJ/kg-K
delH_rxn = 2*10^7; % kJ/kmol
T_feed = 300; % K - It's nat as same as reactor temp
UA = 9000; % kJ/K-s
Tcool = 331; % K

fun = @heat_balance1;

% Stoichiometric matrix
nu = [-1 1];

% Initial value of T
T0 = 310; % K

% Initial value of c
c_0_in = [c0 0]; % kmol/m^3

% Initial Concentration at inlet for feed 
c_0 = c_0_in;
y = [T0 c_0_in];
options = optimoptions('fsolve','Display','iter');
out = fsolve(fun,y,options);

%% post processing 

%conversion of reactant A is : 
% please note there are two "c0" terms :  c_0 and c0, the first one is a
% vector
X = (c0 - out(2))/c0 *100;

T_final = out(1); 

disp("reactor temperature : "+ T_final)
disp("conversion : "+ X)

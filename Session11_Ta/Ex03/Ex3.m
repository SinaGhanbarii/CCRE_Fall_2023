close all; clear variables; clear all;

global b c DiffA DiffB fL HA a KLa KGa kl

%% Data struct
b = 2; %stoichiometric coefficient of B
c = 1; %stoichiometric coefficient of A
DiffA = 1e-6/3600; % diffusion coefficient [m2/s]
DiffB = DiffA; % diffusion coefficient [m2/s]
fL = 0.98; % liquid fraction
HA = 1e5; % [Pa m3/mol]
a = 20; % [m2/m3]
KLa = 20/3600; % mass transfer in the liquid phase [1/s]
KGa = 0.01*(1/3600); % mass transfer in the gas phase [mol/m3/s/Pa]
kl = 1e8/3600; % kinetic constant (liquid volume basis) [m6/mol^2/s]

%% Input data
pTot = 101325; % [Pa]
Ctot = 301; % [mol/m3]
pAin = 5000; % [Pa]
CBin = 1;
CCin = 0;
T = 303; % K
QL = 20/1000; % [m3/s]
QG = 0.8/1000; % [m3/s]
L = 5; % reactor length 
D = 0.40; % diameter [m]

%% Co-current Tower Reactor
FtotG = (pTot/8.314/T)*QG;
FtotL = Ctot*QL;
A = pi/4*D^2;
Vtot = A*L;

%% ODE solution
Yin = [pAin,CBin,CCin]';
[V Y] = ode45(@TowerCoCurrent,[0 Vtot],Yin,[],pTot,Ctot,FtotG,FtotL);

pA = Y(:,1);
CB = Y(:,2);
CC = Y(:,3);
z = V/A;

%% Post-processing
fprintf('A conversion: %f \n', 1-pA(end)/pAin);
figure;
xlabel('axial coordinate [m]');
hold on;
yyaxis left; ylabel('partial pressure of A [Pa]'); 
plot(z,pA);
yyaxis right; ylabel('concentration of B [mol/m3]'); 
plot(z,CB);

figure;
xlabel('axial coordinate [m]');
hold on;
yyaxis left; ylabel('concentration of B [mol/m3]'); 
plot(z,CB);
yyaxis right; ylabel('concentration of C [mol/m3]'); 
plot(z,CC);


%% Reconstructing the overall rate of change and the enhancing factor
for i=1:size(pA,1)
    [r(i),E(i), MH(i), resG(i), resL(i), resK(i)] = OverallRateOfChange(pA(i), CB(i));
end

figure;
plot(z,r*3600);
xlabel('axial coordinate [m]');
ylabel('overall rate of change [mol/m3/hr]'); 


figure;
plot(z,E);
hold on 
plot(z,MH);
xlabel('axial coordinate [m]');
ylabel('enhancing factor'); 
legend('Enhancement factor','Hatta modulus')

figure;
plot(z,resK/3600);
hold on 
plot(z,resG/3600);
hold on 
plot(z,resL/3600);
xlabel('axial coordinate [m]');
ylabel('series resistances'); 
legend('Gas','Film','Reaction')

%% ODE system (defined locally, requires at least MATLAB-2016b)
function dY = TowerCoCurrent(V,Y, pTot, Ctot, FtotG, FtotL)
global b c DiffA DiffB fL HA a KLa KGa kl

pA = Y(1);
CB = Y(2);
r  = OverallRateOfChange(pA,CB);

dpA_over_dV = -pTot/FtotG*r;
dCB_over_dV = -Ctot/FtotL*b*r;
dCC_over_dV = Ctot/FtotL*r;

dY = [dpA_over_dV,dCB_over_dV,dCC_over_dV]';

end

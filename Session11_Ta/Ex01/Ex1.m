clear all; close all; clc;

global b c DiffA DiffB fL HA a KLa KGa kl 

%% Data 
b = 2; %stoichiometric coefficient of B
c = 1; %stoichiometric coefficient of A
DiffA = 1e-6/3600; % diffusion coefficient [m2/s]
DiffB = DiffA; % diffusion coefficient [m2/s]
fL = 0.98; % liquid fraction
HA = 1e5; % [Pa m3/mol]
a = 20; % [m2/m3]
KLa = 20/3600; % mass transfer in the liquid phase [1/s]
KGa = 0.01*(1/3600); % mass transfer in the gas phase [mol/m3/s/Pa]
kl = 1e6/3600; % kinetic constant (liquid volume basis) [m6/mol^2/s]


%% Input data
pAin = 5000; % [Pa]
CBin = 100; %[mol/m3]
pA = pAin; %[Pa]
 
%% Rate computation 
[r,E, Ei, MH, resG, resL, resK] = OverallRateOfChange(pA,CBin);

%% Post processing 
disp("R is : " + r*3600 + " mol/m3/hr")
disp("MH is : " + MH)
disp("E is : " + E)
disp("Ei is : " + Ei)
disp("Series resistances are : " + resG/3600 +", "+ resL/3600 +", "+ resK/3600)

pAinterface = pA - r*resG;
cAinterface = pAinterface/HA;
cAbulk = r*(resK/HA);

disp("Series partial pressures/concentrations of A are : " + pAin +"[Pa], " ... 
    + pAinterface +"[Pa], "+ cAinterface +"[mol/m3], "+ cAbulk +"[mol/m3]")

clear all; close all; clc;

global b c DiffA DiffB fL HA a KLa KGa kl 

%% Data 
b = 1; %stoichiometric coefficient of B
c = 1; %stoichiometric coefficient of A
DiffA = 1e-6; % diffusion coefficient [m2/hr]
DiffB = DiffA; % diffusion coefficient [m2/hr]
fL = 0.1; % liquid fraction
HA = 1; % [Pa m3/mol]
a = 100; % [m2/m3]
KLa = 100; % mass transfer in the liquid phase [1/hr]
KGa = 0.1; % mass transfer in the gas phase [mol/m3/hr/Pa]
kl = 1e+8; % kinetic constant (liquid volume basis) [m3/mol/hr]

pA = 100; %[Pa]
CBin = 100; %[mol/m3]
 
[r,E, MH, resG, resL, resK] = OverallRateOfChange(pA,CBin);
disp("R is : " + r + " mol/m3/hr")
disp("MH is : " + MH)
disp("E is : " + E)
disp("Series resistances are : " + resG +", "+ resL +", "+ resK)

pAinterface = pA - r*resG;
cAinterface = pAinterface/HA;
cAbulk = r*(resK/HA);

disp("Series partial pressures/concentrations of A are : " + pA +"[Pa], " ... 
    + pAinterface +"[Pa], "+ cAinterface +"[mol/m3], "+ cAbulk +"[mol/m3]")


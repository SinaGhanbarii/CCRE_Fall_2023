clear all, clc

%- main script
%
% CCRE - Prof. Matteo Maestri - October 30 - a.a 2023/24
%
%--variables:
%
%--N2 O2 o-xyl - PA - H20 - CO2 - Temp - Press

global Stoichiometry deltaH mw SpecificHeatP Viscosity ...
    TubeDiameter CatalystDensity U ...
    VoidFraction MoltenSaltsTemperature G ParticleDiameter NS NR ...

%% Parameters and Initialization

%Reaction :

Stoichiometry=[0 -3 -1 1 +3 0
               0 -10.5 -1 0 5 8
               0 -7.5 0 -1 2 8];

deltaH=[-1285409 -4564000 -3278591]; %kJ/kmol

%Species:
mw=[28 32 106.16 148.12 18 44]; %kg/kmol
SpecificHeatP=0.992; %kJ/kg/K
Viscosity=2.95e-5; %Pa*s
NS = 6;
NR = 3;

% Reactor:
TubeDiameter=0.0254; %m
ParticleDiameter=0.005; %m
ReactorLength=3; %m
CatalystDensity=2100; %kg/m3
DiameterRatio=TubeDiameter/ParticleDiameter;
VoidFraction=0.363+0.35*(exp(-0.39*DiameterRatio));
U = 385.28; % kJ/m2/h/K

%##################################################################

% Operating Conditions:
MoltenSaltsTemperature=335+273.15; %K
FeedTemperature=MoltenSaltsTemperature;
G=4900; %kg/m2/h
oXylToAirRatio=0.013; 

n_in = [0.79 0.21 oXylToAirRatio 0 0 0]; %moles
ntot = sum(n_in);

for i = 1:NS
    MolarFractionIN(i) = n_in(i)/ntot_in;
end

basis_calc = 1; %kmol
for i= 1:NS
    mass(i) = MolarFractionIN(i)* basis_calc * mw(i);
end
mtot = sum(mass);

for i=1:NS
    MassFractionIN(i) = mass(i)/mtot;
end

FeedPressure = 1.1; %bar


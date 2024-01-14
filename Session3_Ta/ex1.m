%
%--CCRE (Prof. Matteo Maestri) - 2023/2024 - October 9, 2023  
%
% tubular reactor
%
% NO + 0.5 O2 -> NO2  
% Variables: 1) NO || 2) O2 || 3) NO2  || 4) N2 || 5) T
%

clear, clc
global NS nu index DHr_ref Cpmix R Pressure ctot

% Data
NS = 4; % number of species
species = {'NO','O2','NO2','N2'};
index = containers.Map(species,1:NS);
Pressure = 101325; % [Pa]
Temperature = 20+273; % [K]
T_ref = 273; % [K]
Qin_Tref = 1e4/3600; % [m3/s]
Q = Qin_Tref*Temperature/T_ref; % [m3/s]
R = 8.3144; % [J/mol/K]
ctot = Pressure*1e-3/R/Temperature; % [kmol/m3]
xin = [0.09,0.08,0.01,0.82];
mw = [30,32,46,28]; % [kg/kmol]
nu = [-1,-0.5,+1,0]; 
DH = [90200,0,33800,0]; % [KJ/Kg]
Cp = [29.8,29.3,37.9,29.1]; % [KJ/Kg/K]

DHr_ref = 0;

for i = 1:NS
    DHr_ref = DHr_ref + nu(i)*DH(i);
end

Cpmix = 0;
for i = 1:NS
    Cpmix = Cpmix + xin(i)*Cp(i);
end

ntot_in = Q*ctot;
nin = ntot_in*xin;

%% isothermal PFR

[V, nin] = ode23s(@exe1_iso_tubular, 0:0.01:600,nin);

NP = length(V);

for i = 1:NP
    ratio = nin(i, index('NO2'))/nin(i, index('NO'));
    if ratio >= 5
        Vfinal = V(i);
        break
    end
end
disp("Volume of the reactor needed to reach ratio for isothermal case: " + Vfinal + 'm3')

%% Adiabatic Reactor

yin = ones(NS+1,1);
for i=1:NS 
    yin(i) = nin(i);
end

yin(NS+1) = Temperature;
[V, n] = ode23s(@exe1_adi_tubular, 0:0.01:600,yin);
NP = length(V);

for i = 1:NP
    ratio = n(i, index('NO2'))/n(i, index('NO'));
    if ratio >= 5
        Vfinal_adi = V(i);
        break
    end
end
disp("Volume of the reactor needed to reach ratio for adiabatic case: " + Vfinal_adi + 'm3')

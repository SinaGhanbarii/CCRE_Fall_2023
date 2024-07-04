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

%-species:
mw=[28 32 106.16 148.12 18 44]; %kg/kmol
SpecificHeatP=0.992; %kJ/kg/K
Viscosity=2.95e-5; %Pa*s
NS = length(Stoichiometry(1,:));
NR = length(Stoichiometry(:,1));

%-Reactor:
TubeDiameter=0.0254; %m
ParticleDiameter=0.005; %m
ReactorLength=3; %m
CatalystDensity=2100; %kg/m3
DiameterRatio=TubeDiameter/ParticleDiameter;
VoidFraction=0.363+0.35*(exp(-0.39*DiameterRatio));
U = 385.28; % kJ/m2/h/K

%##################################################################

%-Operating Conditions:
MoltenSaltsTemperature=335+273.15; %K
FeedTemperature=MoltenSaltsTemperature;
G=4900; %kg/m2/h
oXylToAirRatio=0.013;
n_in=[0.79 0.21 oXylToAirRatio 0 0 0]; %moles
ntot_in = sum(n_in);
for i = 1:NS
    MolarFractionIN(i) = n_in(i)/ntot_in;
end
FeedPressure=1.3; %bar

%% 
basis_calc=1; %kmol
for i=1:NS
    mass(i)=MolarFractionIN(i)*basis_calc*mw(i); %kg
end
mtot = sum(mass);

for i = 1:NS
    MassFractionIN(i) = mass(i)/mtot;
end

%SYSTEM INTEGRATION:
[t,y]=ode15s('pfr_pseudo_hom',[0:0.01:ReactorLength], ...
    [MassFractionIN FeedTemperature FeedPressure]);

% postprocessing:
NP=length(t);
Across=pi*TubeDiameter^2/4; %m2

for j=1:NP
    for i=1:NS
        n(j,i)=G*Across*y(j,i)/mw(i); % kmol/h
    end
    MolarFraction(j,:)=n(j,:)./sum(n(j,:));
end

%conversion o-xylene:
for j=1:NP
    conv_o_xyl(j)=(n(1,3)-n(j,3))/(n(1,3));
end

%selectivity to PA:
for j=1:NP
    sel_PA(j)=n(j,4)/(n(1,3)-n(j,3)+1e-12);
end

%yield:
for j=1:NP
    yield(j)=conv_o_xyl(j)*sel_PA(j);
end

%productivity: kg_PA/ut in the tube
for j=1:NP
    productivity_tube(j)=yield(j)*n(1,3)*mw(4);
end

target = 8000; % Ton/year
n_tubes = target*1e3/(365*24)*1/productivity_tube(end);
disp("number of tubes: "+ n_tubes)

figure(1)
plot(t,y(:,7)-273.15)
xlabel('Reactor Lenght [m]')
ylabel('Temperature [°C]')
legend('Gas')

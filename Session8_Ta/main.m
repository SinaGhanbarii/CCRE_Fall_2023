clear,clc,clf,close

%- main script
%
% CCRE - Prof. Matteo Maestri - November 13 - a.a 2023/24
%
%--variables:
%
% => N2_bulk O2_bulk o-xyl_bulk - PA_bulk - H20_bulk - CO2_bulk - Temp_bulk 
% => N2_surf O2_surf o-xyl_surf - PA_surf - H20_surf - CO2_surf - Temp_surf 
% => pressure
%
% example: y(2) = y(O2_bulk)
% example: y(2+7) = y(O_2_surf)
%

global Stoichiometry deltaH mw SpecificHeatP Viscosity0 ...
    TubeDiameter CatalystDensity ...
    VoidFraction MoltenSaltsTemperature G ParticleDiameter NS NR lambda0 ...
    lambda_cat ReactorLength TubeThickness lambdaTube ExternalHeatCoefficient

global diffusivity0 a_v_particle NEQ

Stoichiometry=[0 -3 -1 1 +3 0
               0 -10.5 -1 0 5 8
               0 -7.5 0 -1 2 8];

deltaH=[-1285409 -4564000 -3278591]; %kJ/kmol

%-species:
mw=[28 32 106.16 148.12 18 44]; %kg/kmol
SpecificHeatP=0.992; %kJ/kg/K
Viscosity0=2.95e-5; %Pa*s
NS = 6;
NR = 3;

%-Reactor:
TubeDiameter=0.0254; %m
ParticleDiameter=0.005; %m
ReactorLength=3; %m
CatalystDensity=2100; %kg/m3
DiameterRatio=TubeDiameter/ParticleDiameter;
VoidFraction=0.363+0.35*(exp(-0.39*DiameterRatio));

ExternalHeatCoefficient=700; %W/m2/K
TubeThickness=0.0012; %m
lambdaTube=20; %W/m/K 
lambda0=4.78e-2; % W/m/K
lambda_cat=1.5; %W/m/K

diffusivity0 = [8.31e-05, 7.02e-05, 2.58e-05, 2.22e-05, 8.98e-05, 2.02e-05];
a_v_particle = 6/ParticleDiameter*(1-VoidFraction);

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

basis_calc=1; %kmol
for i=1:NS
    mass(i)=MolarFractionIN(i)*basis_calc*mw(i); %kg
end
mtot = sum(mass);

for i = 1:NS
    MassFractionIN(i) = mass(i)/mtot;
end

NEQ = 2*NS+2+1; 
M = zeros(NEQ,NEQ);
for i = 1:(NS+1)
    M(i,i) = 1;
end
M(NEQ,NEQ) = 1; % pressure

options = odeset('Mass',M,'AbsTol',1e-06,'RelTol',1e-06);

[t,y]=ode15s(@pfr_het,[0:0.001:ReactorLength], ...
    [MassFractionIN FeedTemperature ...
    MassFractionIN FeedTemperature FeedPressure],options);

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
hold on
plot(t,y(:,end-1)-273.15,'Linestyle','--')
xlabel('Reactor Lenght [m]')
ylabel('Temperature [°C]')
legend('Gas','Surface')

figure(2)
plot(t,y(:,3))
hold on
plot(t,y(:,3+7),'linestyle','--')
xlabel('Reactor Lenght [m]')
ylabel('O-xylene mass fraction [-]')
legend('Gas','Surface')

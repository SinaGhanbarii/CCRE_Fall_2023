clear all, clc

%- main script
%
% Procedure: 
%   1: Declare how to answer the question!
%   2: Explain the procedure that will lead you to answer the question 
%   3: Declare and justify the assumptions that you are going to make 
%   4: Report equations/models ---> Write Boundary conditions and explain
%   each parameter and physical meaning of each term!
%   5: Report all the transformations in terms of the units of measurement in order to have all the parameters 
%   coherent
%   6: report the results. If it is a graph, report a 
%   qualitative curve with quantitative values of the quantity at given points (e.g., maximum, minimum, 
%   outlet, …). 
% CCRE - Prof. Matteo Maestri - November 13 - a.a 2023/24
%
%--variables:
%
% => N2_bulk O2_bulk o-xyl_bulk - PA_bulk - H20_bulk - CO2_bulk - Temp_bulk 
% => N2_surf O2_surf o-xyl_surf - PA_surf - H20_surf - CO2_surf - Temp_surf 
% => pressure
% CH3OH O2 CH2O H2O CO N2 T P 
% 
% example: y(2) = y(O2_bulk)
% example: y(2+7) = y(O_2_surf)
%

global Stoichiometry deltaH mw SpecificHeatP Viscosity0 ...
    TubeDiameter CatalystDensity ...
    VoidFraction MoltenSaltsTemperature G ParticleDiameter NS NR lambda0 ...
    lambda_cat ReactorLength TubeThickness lambdaTube ExternalHeatCoefficient ...
    Sphericity

global diffusivity0 a_v_particle NEQ

Stoichiometry=[-1  -0.5  +1   +1  0  0
                0  -0.5  -1   +1  +1 0];

deltaH=[-38 -57]*4186*1000; %kJ/kmol

%-species:
mw=[32 32 30 18 28 28]; %kg/kmol
SpecificHeatP=9.5*4186; %kJ/kmol/K
Viscosity0=5.4e-2/3600; %Pa*s
NS = 6;
NR = 2;

%-Reactor:
TubeDiameter=0.016; %m
ParticleDiameter=0.003; %m
Sphericity = 0.8;
ReactorLength=0.8; %m
CatalystDensity=1450; %kg/m3
DiameterRatio=TubeDiameter/ParticleDiameter;
VoidFraction=0.55;

ExternalHeatCoefficient=600*4186/3600; %W/m2/K
TubeThickness=0.0015; %m
lambdaTube=18*4186/3600; %W/m/K 
lambda0=(2.3e-2)*4186/3600; % W/m/K
lambda_cat=1.318*4186/3600; %W/m/K

diffusivity0 = [0.1 0.1 0.1 0.1 0.1 0.1]/3600; %m2/s
a_v_particle = 6/ParticleDiameter*(1-VoidFraction);

%-Operating Conditions:
MoltenSaltsTemperature=285+273.15; %K
FeedTemperature=245+273.15; %K
% G=4900; %kg/m2/h
MolarFlowratePerTube = 1200*101.324/(8.314*273.15)/3600/1000; %kmol/s
MetanolToAirRatio=0.085;
n_in=[MetanolToAirRatio 0.21 0 0 0 0.79]; %moles
InletMolarFlowRate = MolarFlowratePerTube.*n_in; %kmol/s
InletMassFlowRate = InletMolarFlowRate.*mw; %kg/s
TotalMassFlowRate = sum(InletMassFlowRate); %kg/s 
MassFractionIN = InletMassFlowRate./TotalMassFlowRate;

Across=pi*TubeDiameter^2/4; %m2
G = (TotalMassFlowRate/Across); %kg/s/m2

% ntot_in = sum(n_in);
% for i = 1:NS
%     MolarFractionIN(i) = n_in(i)/ntot_in;
% end
FeedPressure=1.2*1.01325; %bar
% 
% basis_calc=1; %kmol
% for i=1:NS
%     mass(i)=MolarFractionIN(i)*basis_calc*mw(i); %kg
% end
% mtot = sum(mass);
% 
% for i = 1:NS
%     MassFractionIN(i) = mass(i)/mtot;
% end

% Mass matrix
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

% target = 8000; % Ton/year
% n_tubes = target*1e3/(365*24)*1/productivity_tube(end);
% disp("number of tubes: "+ n_tubes)

Conversion = (y(1,1)-y(end,1))/y(1,1)   % X
Selectivity = -1/1 * (y(3,end)-y(3,1))/(y(1,1)-y(1,end))    % Sigma
Yield = Conversion*Selectivity          % Zeta

figure(1)
plot(t,y(:,7)-273.15)
hold on
plot(t,y(:,end-1)-273.15,'Linestyle','--')
xlabel('Reactor Lenght [m]')
ylabel('Temperature [°C]')
legend('Gas','Surface')

figure(2)
plot(t,y(:,1))
hold on
plot(t,y(:,1+7),'linestyle','--')
xlabel('Reactor Lenght [m]')
ylabel('O-xylene mass fraction [-]')
legend('Gas','Surface')


% Sherwood Number (Sh): Ratio of convective mass transfer to diffusive mass transfer.
% Schmidt Number (Sc): Ratio of momentum diffusivity (kinematic viscosity) to mass diffusivity.
% Reynolds Number (Re): Ratio of inertial forces to viscous forces.
% Nusselt Number (Nu): Ratio of convective to conductive heat transfer.
% Prandtl Number (Pr): Ratio of momentum diffusivity (kinematic viscosity) to thermal diffusivity.
% Peclet Number (Pe): Ratio of convective to diffusive transport rates, used in both heat and mass transfer.

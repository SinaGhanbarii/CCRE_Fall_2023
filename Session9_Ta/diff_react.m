clear all, clc

%- main script
%
% CCRE - Prof. Matteo Maestri - November 20 - a.a 2023/24
%
%--variables:
%
%--N2 O2 o-xyl - PA - H20 - CO2 - Temp

global NR NS NEQ mw N_reactions deltaH Stoichiometry diffusivity0 epsi tau pore_diameter ...
    delta_r CatalystDensity pressure r_vector lambda_cat yIN


% i => index for the species
% j => index for the radial points

%% Parameters and Initialization
%Reaction:
Stoichiometry=[0 -3 -1 1 +3 0
               0 -10.5 -1 0 5 8
               0 -7.5 0 -1 2 8];

deltaH=[-1285409 -4564000 -3278591]; %kJ/kmol
diffusivity0=[8.31e-5 7.02e-5 2.58e-5 2.22e-5 8.98e-5 2.02e-5]; %m2/s

%-species:
mw=[28 32 106.16 148.12 18 44]; %kg/kmol
NS = 6; % number of species
N_reactions = 3; % number of reactions
NEQ = NS+1; % number of equations

%-Reactor:
ParticleDiameter=0.005; %m
ParticleRadius = ParticleDiameter/2;
CatalystDensity=2100; %kg/m3
lambda_cat=1.5; %W/m/K

Temperature=335+273.15; %K
oXylToAirRatio=0.013;
n_in=[0.79 0.21 oXylToAirRatio 0 0 0]; %moles
ntot_in = sum(n_in);
for i = 1:NS
    MolarFractionIN(i) = n_in(i)/ntot_in;
end
pressure=1.3*1e5; %Pa

basis_calc=1; %kmol
for i=1:NS
    mass(i)=MolarFractionIN(i)*basis_calc*mw(i); %kg
end
mtot = sum(mass);

for i = 1:NS
    MassFractionIN(i) = mass(i)/mtot;
end

%-PELLET PARAMETERS:
NR=80; % Number of radial points
initial_point = 0e-3;
r_vector=linspace(initial_point,ParticleRadius,NR);
delta_r = r_vector(2)-r_vector(1);

epsi=0.3;
tau=5;
pore_diameter=1e-8; %m

%initial guess for the fsolve
yIN = [MassFractionIN Temperature];
y_trial = zeros(1,NR*NEQ); 

for i = 1:NEQ
    for j = 1:NR
        y_trial((j-1)*NEQ+i) = yIN(i);
    end
end

y = fsolve('diffusion_reaction',y_trial); 

%% postprocessing:
solution = zeros(NR,NEQ);

for i = 1:NEQ
    for j = 1:NR
        solution(j,i) = y((j-1)*NEQ+i); 
    end
end

mol=zeros(NR,NS);
basis_calc = 1; % kgtot
for j=1:NR
    for i=1:NS
        mol(j,i)=solution(j,i)/mw(i)*basis_calc; %mol
    end
    MolarFractionSolid(j,:)=mol(j,:)./sum(mol(j,:));
end

%conversion o-xylene:
for j=1:NR
    conv_o_xyl(j)=(mol(NR,3)-mol(j,3))/(mol(NR,3));
end

%selectivity to PA:
for j=1:NR
    sel_PA(j)=mol(j,4)/(mol(NR,3)-mol(j,3)+1e-12);
end

%yield:
for j=1:NR
    yield(j)=conv_o_xyl(j)*sel_PA(j);
end

% evaluation of the catalyst efficiency
TemperatureSolid = solution(:,end);
Pressure_bar = pressure/1e5; % bar
NetRateProduction=zeros(NR,NS);
RateH=zeros(NR,1);

for j=1:NR
    ReactionK(1)=exp(19.837-13636/TemperatureSolid(j));
    ReactionK(2)=exp(18.970-14394/TemperatureSolid(j));
    ReactionK(3)=exp(20.860-15803/TemperatureSolid(j));
    ReactionRate(1)=ReactionK(1)*Pressure_bar*MolarFractionSolid(j,3)*...
        Pressure_bar*MolarFractionSolid(j,2);
    ReactionRate(2)=ReactionK(2)*Pressure_bar*MolarFractionSolid(j,3)*...
        Pressure_bar*MolarFractionSolid(j,2);
    ReactionRate(3)=ReactionK(3)*Pressure_bar*MolarFractionSolid(j,4)*...
        Pressure_bar*MolarFractionSolid(j,2); %kmol/kg_cat_/h

    for k=1:N_reactions                    
        for i=1:NS
            NetRateProduction(j,i)=NetRateProduction(j,i)+...
                Stoichiometry(k,i)*ReactionRate(k);
        end
        RateH(j)=RateH(j)+deltaH(k)*ReactionRate(k); %kJ/kg_cat/h
    end
end

r_vector = (r_vector-initial_point)/(ParticleRadius-initial_point);

V=4/3*pi*r_vector.^3;
Ri=trapz(V,NetRateProduction(:,3)*CatalystDensity*mw(3)/3600)/(4/3*pi*(r_vector(end)^3-r_vector(1)^3));
Rsurface = NetRateProduction(end,3)*CatalystDensity*mw(3)/3600;

catalyst_efficiency = Ri/Rsurface*100;
disp("eta: "+catalyst_efficiency)

Ri_actual = zeros(NR,1);
for j = 1:NR
    Ri_actual(j) = NetRateProduction(j,3)*CatalystDensity*mw(3)/3600;
end

effectiveness_factor = Ri_actual./Rsurface;

% ------
figure(1)
plot(r_vector,solution(:,7)-273.15)
hold on
xlabel('Pellet radius [m]')
ylabel('Temperature [°C]')

figure(2)
plot(r_vector,solution(:,3))
hold on
xlabel('Pellet radius [m]')
ylabel('O-xylene [-]')

figure(3)
plot(r_vector,effectiveness_factor)
hold on
xlabel('Pellet radius [m]')
ylabel('catalyst efficiency [-]')

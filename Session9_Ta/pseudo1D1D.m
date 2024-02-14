clear all

%
%-- CCRE - Prof. Matteo Maestri - November 20 - a.a 2023/24  
%
% => N2_bulk O2_bulk o-xyl_bulk PA_bulk H2O_bulk CO2_bulk Temperature_bulk 
% => N2_solid O2_solid o-xyl_solid PA_solid H2O_solid CO2_solid Temperature_solid
%
% => N2_pellet_radial_profile O2_pellet_radial_profile o-xyl_pellet_radial_profile 
% => PA_pellet_radial_profile H2O_pellet_radial_profile CO2_pellet_radial_profile 
% => Temperature_pellet_radial_profile
% => Pressure
%
% example: y(2)=omega(O2_bulk)
% example: y(2+7)=omega(O2_solid)
% example: y(startSolidIndex+(j-1)*NEQ+2)=omega(O2_jth_radial_position)
%

global Stoichiometry deltaH TubeDiameter ParticleDiameter ReactorLength...
    CatalystDensity VoidFraction MoltenSaltsTemperature...
    G SpecificHeatP Viscosity0 mw

global ExternalHeatCoefficient TubeThickness lambdaTube lambda0 lambda_cat
global diffusivity0 a_v_particle
global NS NEQ NR r_vector startSolidIndex epsi tau pore_diameter
global initial_point ParticleRadius

%% PROBLEM DATA
%Reactions: 

Stoichiometry=[0 -3 -1 +1 +3 0
               0 -10.5 -1 0 5 8
               0 -7.5 0 -1 2 8];
           
deltaH=[-1285409 -4564000 -3278591]; %kJ/mol

%Species:
mw=[28 32 106.16 148.12 18 44]; % kg/kmol

ParticleDiameter=0.005; %m
ParticleRadius = ParticleDiameter/2;
NR=50; % Number of radial points
initial_point = 0e-3;
r_vector=linspace(initial_point,ParticleRadius,NR);
delta_r = r_vector(2)-r_vector(1);

epsi=0.3;
tau=5;
pore_diameter=1e-8; %m

%-REACTOR PARAMETERS:
TubeDiameter=0.0254; %m 
ReactorLength=3; %m 
CatalystDensity=2100; %kg/m3
DiameterRatio=TubeDiameter/ParticleDiameter;
VoidFraction=0.363+0.35*(exp(-0.39*DiameterRatio)); 

ExternalHeatCoefficient=700; %W/m2/K
TubeThickness=0.0012; %m
lambdaTube=20; %W/m/K 
lambda0=4.78e-2; % W/m/K
lambda_cat=1.5; %W/m/K

diffusivity0=[8.31e-5 7.02e-5 2.58e-5 2.22e-5 8.98e-5 2.02e-5]; %m2/s
a_v_particle=6*(1-VoidFraction)/ParticleDiameter;

%-OPERATING CONDITIONS:
MoltenSaltsTemperature=335+273.15; %K %360 %380
FeedTemperature=MoltenSaltsTemperature; 
G=4900; %kg/m2/h
oXylToAirRatio=0.013; %0.013
n_in=[0.79 0.21 oXylToAirRatio 0 0 0]; 
MolarFractionIN=n_in./sum(n_in);
SpecificHeatP=0.992; %kJ/kg/K
Viscosity0=2.95e-5; %Pa*s
FeedPressure=1.3; %bar
NS = 6; % number of species

%==
basis_calc=1; %kmol YOUR CHOICE!
for i=1:NS
    mass(i)=basis_calc*MolarFractionIN(i)*mw(i); %kg
end
MassFractionIN=mass./sum(mass);

%% Generation of y0 vector
NEQ=NS+1;
y0=zeros((NR+1)*NEQ+1,1); % y0=zeros(NS*NR+NR+NS+1+1,1);
startSolidIndex=NEQ;
for i=1:NS
   y0(i)=MassFractionIN(i); % initial conditions
   for j=1:NR
       y0(startSolidIndex+(j-1)*NEQ+i)=MassFractionIN(i); % initial guess
   end
end

i=NEQ;
y0(i)=FeedTemperature; % initial conditions
for j=1:NR
    y0(startSolidIndex+(j-1)*NEQ+i)=FeedTemperature; % initial guess
end

y0(end)=FeedPressure;

%% Heterogeneous Mass Matrix
%solution of the ODE-1st-order-system
%===============================================
%-- Mass Matrix: it defines the DAE system
%
M=zeros((NR+1)*NEQ+1,(NR+1)*NEQ+1);
for i=1:NS+1
    M(i,i)=1;
end
M(end,end)=1;
%
options=odeset('Mass',M,'AbsTol',1e-6,'RelTol',1e-6);

%% SYSTEM INTEGRATION
[t,y]=ode15s(@PFR_1D1D_pseudo,...
  [0:0.001:ReactorLength],y0,options);

%% POST PROCESSING:

NP=size(y(:,1)); %number of points
NP=NP(1);

Across=pi*TubeDiameter^2/4;

for j=1:NP
    for i=1:NS
        n(j,i)=G*Across*y(j,i)/mw(i); %kmol/h
    end
    MolarFraction(j,:)=n(j,:)./sum(n(j,:));
end

for j=1:NP
    conv_oxyl(j)=(n(1,3)-n(j,3))/n(1,3);
end

for j=1:NP
    sel_PA(j)=n(j,4)/(n(1,3)-n(j,3)+1d-12);
end

for j=1:NP
    yield(j)=sel_PA(j)*conv_oxyl(j);
end

for j=1:NP
    prod_per_tube(j)=yield(j)*n(1,3)*mw(4); %kg_PA/h/tube
end

oxy_1D_1D=zeros(NP,NR);
T_1D_1D=zeros(NP,NR);
co2_1D_1D=zeros(NP,NR);

for k=1:NP 
    for j=1:NR
        oxy_1D_1D(k,j)=y(k,startSolidIndex+(j-1)*NEQ+3);
        co2_1D_1D(k,j)=y(k,startSolidIndex+(j-1)*NEQ+6);
        T_1D_1D(k,j)=y(k,startSolidIndex+(j-1)*NEQ+NEQ);
    end
end

%% 1D PLOT
figure(1)
plot(t,y(:,7)-273.15)
hold on
plot(t,y(:,end-1)-273.15,'linestyle','--')
xlabel('Reactor Length [m]')
ylabel('Temperature [°C]')
legend('Gas','Surface')

figure(2)
plot(t,y(:,end))
xlabel('Reactor Length [m]')
ylabel('Pressure [bar]')

figure(3)
jj=NR;
plot(t,y(:,3))
hold on
plot(t,y(:,startSolidIndex+(jj-1)*NEQ+3),'linestyle','--')
xlabel('Reactor Length [m]')
ylabel('O-xylene mass fraction [-]')
legend('Gas','Surface')

%% 2D plot
[Z,R]=meshgrid(t,r_vector);

figure(4)
s=surf(Z,R,oxy_1D_1D');
colormap('jet')
s.EdgeColor = 'none';
xlabel('Reactor Length [m]')
ylabel('Pellet Radius')
zlabel('Catalyst o-xylene mass fraction[-]')
xlim([0 ReactorLength])
ylim([0 ParticleDiameter/2])    

figure(5)
s=surf(Z,R,T_1D_1D'-273.15);
colormap('jet')
s.EdgeColor = 'none';
xlabel('Reactor Length [m]')
ylabel('Pellet Radius')
zlabel('Catalyst Temperature [°C]')
xlim([0 ReactorLength])
ylim([0 ParticleDiameter/2])    

figure(6)
s=surf(Z,R,co2_1D_1D');
colormap('jet')
s.EdgeColor = 'none';
xlabel('Reactor Length [m]')
ylabel('Pellet Radius')
zlabel('Catalyst {CO}_2 mass fraction[-]')
xlim([0 ReactorLength])
ylim([0 ParticleDiameter/2])    

result=sprintf('Productivity per tube is: %.3d [kg/h]',prod_per_tube(end));
disp(result);

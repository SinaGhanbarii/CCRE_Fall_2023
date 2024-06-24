function pfr_het=pfr_het(t,y)

%
% CCRE - Prof. Matteo Maestri - November 13 - a.a 2023/24
%

global Stoichiometry deltaH mw SpecificHeatP Viscosity0 ...
    TubeDiameter CatalystDensity ...
    VoidFraction MoltenSaltsTemperature G ParticleDiameter NS NR lambda0 ...
    lambda_cat ReactorLength TubeThickness lambdaTube ExternalHeatCoefficient ...
    Sphericity

global diffusivity0 a_v_particle NEQ

Nsolid = NS+1;

for i=1:NS
    mol(i)=y(i)/mw(i); %mol
    molSurface(i)=y(i+Nsolid)/mw(i);
end
MolarFraction=mol./sum(mol);
MolarFractionSurface=molSurface./sum(molSurface);

Temperature=y(Nsolid);
TemperatureSurface=y(end-1);
Pressure_bar=y(end);
Pressure_Pa=Pressure_bar*1e5;
Pressure_atm = Pressure_bar/1.01325;
DensityMolarGas=1e-3*Pressure_Pa/8.314/Temperature; %kmol/m3

mw_av=0;
for i=1:NS
    mw_av=mw_av+MolarFraction(i)*mw(i);
end

DensityMassGas=DensityMolarGas*mw_av; %kg/m3

SuperficialVelocity=G/DensityMassGas; %m/s
SpecificHeatkg=SpecificHeatP/mw_av;

%-KINETIC SCHEME: CH3OH O2 CH2O H2O CO N2 T P 
% -- 

PartialPressure = Pressure_atm*MolarFractionSurface;
pCH3OH = PartialPressure(1);
pO2 = PartialPressure(2);
pCH2O = PartialPressure(3);

eta1= 0.82*((TemperatureSurface-273.15)/230)^(-3.3);
eta2 = 0.75*((TemperatureSurface-273.15)/230)^(-3);

eta = [eta1 eta2];
disp(eta)
R=8.314*1000/4186;  %kcal/kmol

ReactionK(1)=exp(18.15-17700/(TemperatureSurface*R));
ReactionK(2)=2.42e5*exp(-16050/(TemperatureSurface*R));
ReactionRate(1)=ReactionK(1)*(sqrt(pO2*pCH3OH)/(pO2^0.5+0.5*pCH3OH^0.5))*101.325/8.314/273.15*1000/3600/1000; %kmol/kg_cat_/h
ReactionRate(2)=ReactionK(2)*pCH2O*101.325/8.314/273.15*1000/3600/1000; %kmol/kg_cat_/h

for i=1:NS
NetRateProduction(i)=0;
end

for j=1:NR
    for i=1:NS
        NetRateProduction(i)=NetRateProduction(i)+...
                              Stoichiometry(j,i)*ReactionRate(j)*eta(j);
    end
end

%% EXTERNAL HEAT TRANSFER
lambda=lambda0;
Viscosity=Viscosity0;

Re=G*ParticleDiameter/Viscosity; 
Prandtl=Viscosity*SpecificHeatkg/lambda; 

%--Dixon + Specchia:-----
lambda_static=lambda*(VoidFraction+...
     (1-VoidFraction)/(0.22*VoidFraction^2+2/3*(lambda/lambda_cat)));

Pe_rf=8.65*(1+19.4*(ParticleDiameter/TubeDiameter)^2);
lambda_dynamic=lambda*(Re*Prandtl/Pe_rf);

lambda_eff=lambda_static+lambda_dynamic;

Static=2*VoidFraction...
 +(1-VoidFraction)/(0.0024*(TubeDiameter/ParticleDiameter)^1.58...
 +(1/3)*(lambda/lambda_cat));
alfa_w_static=Static*lambda/ParticleDiameter; %W/m2/K
if Re<1200
    alfa_w_dynamic=(lambda/ParticleDiameter)*(0.0835*Re^0.91);
end
if Re>=1200
    alfa_w_dynamic=(lambda/ParticleDiameter)*(1.23*Re^0.53);
end
alfa_w=alfa_w_dynamic+alfa_w_static;

Bi=alfa_w*TubeDiameter/lambda_eff;
A=6*(Bi+4)/(Bi+3);
InternalHeatCoefficient=alfa_w/(1+Bi/A); %W/m2/K

Ai=pi*TubeDiameter*ReactorLength;
Ae=pi*(TubeDiameter+2*TubeThickness)*ReactorLength;
Aln=(Ae-Ai)/log(Ae/Ai);

U=(1/InternalHeatCoefficient+TubeThickness/lambdaTube*Ai/Aln+1/ExternalHeatCoefficient*Ai/Ae)^-1; %W/m2/K; 

%% INTERPHASE Trasport
Re_p = (1/Sphericity)*(G)*ParticleDiameter/6/Viscosity/(1-VoidFraction);
J_m  = Sphericity*0.61*Re_p^(-0.41); % YOSHIDA CORRELATION
J_h = J_m; % Chilton-Colburn Analogy 

Re = (G)*ParticleDiameter/Viscosity;
Nusselt = J_h*Re*Prandtl^(1/3);

h_thermal = Nusselt*lambda0/ParticleDiameter;

for i = 1:NS
    Sc(i) = Viscosity/DensityMassGas/diffusivity0(i);
    Sh(i) = J_m*Re*Sc(i)^(1/3);
    k_mat(i) = Sh(i)*diffusivity0(i)/ParticleDiameter;
end

%% -- Governing Equations 
for i=1:NS
    pfr_het(i)=(a_v_particle*k_mat(i)*DensityMassGas*(y(i+Nsolid)-y(i)))/G;
    pfr_het(i+Nsolid)=-(a_v_particle*k_mat(i)*DensityMassGas*(y(i+Nsolid)-y(i)))+ ...
        NetRateProduction(i)*(1-VoidFraction)*CatalystDensity*mw(i);
end

RateH=0;
for j=1:NR
    RateH=RateH+deltaH(j)*ReactionRate(j)*eta(j); %kJ/kg_cat/h
end

pfr_het(Nsolid)=(h_thermal*a_v_particle*(TemperatureSurface-Temperature)+ ...
    U*(4/TubeDiameter)*(MoltenSaltsTemperature-Temperature))/G/SpecificHeatkg;

pfr_het(NEQ-1)=(-RateH*(1-VoidFraction)*CatalystDensity - ...
    h_thermal*a_v_particle*(TemperatureSurface-Temperature))/G/SpecificHeatkg;

pfr_het(NEQ)=-(((150*(1-VoidFraction)^2/VoidFraction^3))*...
    Viscosity*SuperficialVelocity/ParticleDiameter^2+...
    7/4*(DensityMassGas*SuperficialVelocity^2/ParticleDiameter)...
    *((1-VoidFraction)/VoidFraction^3))/1e5;

pfr_het=pfr_het';
end

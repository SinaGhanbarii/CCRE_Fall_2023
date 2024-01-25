function pfr_pseudo_hom=pfr_pseudo_hom(t,y)

%
% CCRE - Prof. Matteo Maestri - November 8 - a.a 2023/24
%

global Stoichiometry deltaH mw SpecificHeatP Viscosity0 ...
    TubeDiameter CatalystDensity ...
    VoidFraction MoltenSaltsTemperature G ParticleDiameter NS NR lambda0 ...
    lambda_cat ReactorLength TubeThickness lambdaTube ExternalHeatCoefficient


basis_calc = 1; % kgtot

for i=1:NS
    mol(i)=y(i)/mw(i)*basis_calc; %mol
end
MolarFraction=mol./sum(mol);

Temperature=y(end-1);
Pressure_bar=y(end);
Pressure_Pa=Pressure_bar*1e5;
DensityMolarGas=1e-3*Pressure_Pa/8.314/Temperature; %kmol/m3

mw_av=0;
for i=1:NS
    mw_av=mw_av+MolarFraction(i)*mw(i);
end

DensityMassGas=DensityMolarGas*mw_av; %kg/m3

SuperficialVelocity=G/DensityMassGas/3600; %m/s

%-KINETIC SCHEME:
ReactionK(1)=exp(19.837-13636/Temperature);
ReactionK(2)=exp(18.970-14394/Temperature);
ReactionK(3)=exp(20.860-15803/Temperature);
ReactionRate(1)=ReactionK(1)*Pressure_bar^2....
    *MolarFraction(2)*MolarFraction(3); %kmol/kg_cat_/h
ReactionRate(2)=ReactionK(2)*Pressure_bar^2....
    *MolarFraction(2)*MolarFraction(3); %kmol/kg_cat_/h
ReactionRate(3)=ReactionK(3)*Pressure_bar^2....
    *MolarFraction(2)*MolarFraction(4); %kmol/kg_cat_/h

%% HEAT TRANSFER
lambda=lambda0;
Viscosity=Viscosity0;

Re=G/3600*ParticleDiameter/Viscosity; 
Prandtl=Viscosity*SpecificHeatP/lambda*1e3; 

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

U=(1/InternalHeatCoefficient+TubeThickness/lambdaTube*Ai/Aln+1/ExternalHeatCoefficient*Ai/Ae)^-1; %W/m2/K
U = U*3600/1000; % kJ/m2/h/K


%% -- Governing Equations (TO DO)
NetRateProduction = zeros(NS,1);

for j=1:NR
    for i=1:NS
        NetRateProduction(i)=NetRateProduction(i)+...
                              Stoichiometry(j,i)*ReactionRate(j);
    end
end

% mass balances
for i=1:NS
    pfr_pseudo_hom(i)=...
        NetRateProduction(i)*(1-VoidFraction)*CatalystDensity*mw(i)/G;
end

RateH=0;
for j=1:NR
    RateH=RateH+deltaH(j)*ReactionRate(j); %kJ/kg_cat/h
end

% energy balance
pfr_pseudo_hom(NS+1)=(-RateH*(1-VoidFraction)*CatalystDensity+...
    U*(4/TubeDiameter)...
    *(MoltenSaltsTemperature-Temperature))/G/SpecificHeatP;

% ergun equation
pfr_pseudo_hom(NS+2)=-(((150*(1-VoidFraction)^2/VoidFraction^3))*...
    Viscosity*SuperficialVelocity/ParticleDiameter^2+...
    7/4*(DensityMassGas*SuperficialVelocity^2/ParticleDiameter)...
    *((1-VoidFraction)/VoidFraction^3))/1e5;

pfr_pseudo_hom=pfr_pseudo_hom';
end


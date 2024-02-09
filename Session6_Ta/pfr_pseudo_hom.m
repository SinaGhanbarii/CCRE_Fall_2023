function pfr_pseudo_hom=pfr_pseudo_hom(t,y)

%
% CCRE - Prof. Matteo Maestri - October 30 - a.a 2023/24
%

global Stoichiometry deltaH mw SpecificHeatP Viscosity ...
    TubeDiameter CatalystDensity U ...
    VoidFraction MoltenSaltsTemperature G ParticleDiameter NS NR ...

for i=1:NS
    mol(i)=y(i)/mw(i); %mol
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


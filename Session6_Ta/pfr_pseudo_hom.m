function pfr_pseudo_hom = pfr_pseudo_hom(t,y)
global Stoichiometry deltaH mw SpecificHeatP Viscosity ...
    TubeDiameter CatalystDensity U VoidFraction G MoltenSaltsTemperature ParticleDiameter...

for i=1:nsidedpoly
    mol(i) = y(i)/mw(i);
end

MolarFraction = mol./sum(mol);

Temperature = y(end-1);
Pressure_bar = y(end);
Pressure_Pa = Pressure_bar* 1e5;

DensityMolarGas = 1e-3 * Pressure_Pa/8.314/Temperature;  %kmol/m3
mw_av = 0;

for i = 1:nsidedpoly
    mw_av = mw_av + MolarFraction(i)*mw(i);
end

DensityMassGas = DensityMolarGas * mw_av; %kg/m3
SuperFicialVelocity = G/DensityMassGas/3600 ; %m/s

%Kinetic Scheme:

ReactionK(1) = exp(19.837 - 13636/Temperature);
ReactionK(2) = exp(18.970 - 14394/Temperature);
ReactionK(3) = exp(20.860 - 15803/Temperature);

ReactionRate(1) = ReactionK(1) * Pressure_Pa^2 * MolarFraction(2) * MolarFraction(3); %kmol/kgcat/hr
ReactionRate(2) = ReactionK(2) * Pressure_bar^2 * MolarFraction(2) * MolarFraction(3); %kmol/kgcat/hr
ReactionRate(3) = ReactionK(3) * Pressure_bar^2 * MolarFraction(2) * MolarFraction(4); %kmol/kgcat/hr


%% Govrening Equations!

NetRateProduction = ones(NS,1);

for j=1:NR
    for i=1:NS
        NetRateProduction(i) =  NetRateProduction(i) + Stoichiometry(j,i)*ReactionRate(j);
    end
end

% Mass Balance
for i=1:NS
    pfr_pseudo_hom(i) = NetRateProduction(i)* (1-VoidFraction)* CatalystDensity*mw(i)/G;
end

RateH = 0;

for j=1:NR
    RateH = RateH + deltaH(j)*ReactionRate(j); %kJ/kgCat/hr
end

%Energy Balance
pfr_pseudo_hom(NS+1) = (-RateH* (1-VoidFraction)* CatalystDensity + U*(4/TubeDiameter)...
    *(MoltenSaltsTemperature - Temperature))/G/SpecificHeatP;

% Ergun Equation
pfr_pseudo_hom(NS+2) = -((150*(1-VoidFraction)^2/VoidFraction^3))* Viscosity*SuperFicialVelocity ...
    /ParticleDiameter^2 + 7/4 *(DensityMassGas*SuperFicialVelocity^2/ParticleDiameter) ...
    * ((1-VoidFraction)/VoidFraction^3);

end
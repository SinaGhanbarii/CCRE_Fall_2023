function pfr_pseudo_hom = pfr_pseudo_hom(t,y)
global Stoichiometry deltaH mw SpecificHeatP Viscosity ...
    TubeDiameter CatalystDensity U ...

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
SuperFicialVelocity = G/DensityMassGas;




end
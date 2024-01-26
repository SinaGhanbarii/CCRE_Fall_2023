function pfr_het=pfr_het(t,y)

%
% CCRE - Prof. Matteo Maestri - November 13 - a.a 2023/24
%

global Stoichiometry deltaH mw SpecificHeatP Viscosity0 ...
    TubeDiameter CatalystDensity ...
    VoidFraction MoltenSaltsTemperature G ParticleDiameter NS NR lambda0 ...
    lambda_cat ReactorLength TubeThickness lambdaTube ExternalHeatCoefficient

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
DensityMolarGas=1e-3*Pressure_Pa/8.314/Temperature; %kmol/m3

mw_av=0;
for i=1:NS
    mw_av=mw_av+MolarFraction(i)*mw(i);
end

DensityMassGas=DensityMolarGas*mw_av; %kg/m3

SuperficialVelocity=G/DensityMassGas/3600; %m/s

%-KINETIC SCHEME:
% -- 
ReactionK(1)=exp(19.837-13636/TemperatureSurface);
ReactionK(2)=exp(18.970-14394/TemperatureSurface);
ReactionK(3)=exp(20.860-15803/TemperatureSurface);
ReactionRate(1)=ReactionK(1)*Pressure_bar^2....
    *MolarFractionSurface(2)*MolarFractionSurface(3); %kmol/kg_cat_/h
ReactionRate(2)=ReactionK(2)*Pressure_bar^2....
    *MolarFractionSurface(2)*MolarFractionSurface(3); %kmol/kg_cat_/h
ReactionRate(3)=ReactionK(3)*Pressure_bar^2....
    *MolarFractionSurface(2)*MolarFractionSurface(4); %kmol/kg_cat_/h

for i=1:NS
NetRateProduction(i)=0;
end

for j=1:NR
    for i=1:NS
        NetRateProduction(i)=NetRateProduction(i)+...
                              Stoichiometry(j,i)*ReactionRate(j);
    end
end

%% EXTERNAL HEAT TRANSFER
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

%% INTERPHASE Trasport
Re_p = (G/3600)*ParticleDiameter/6/Viscosity/(1-VoidFraction);
J_m  =0.61*Re_p^(-0.41); % YOSHIDA CORRELATION
J_h = J_m; % Chilton-Colburn Analogy 

Re = (G/3600)*ParticleDiameter/Viscosity;
Nusselt = J_h*Re*Prandtl^(1/3);

h_thermal = Nusselt*lambda0/ParticleDiameter;

for i = 1:NS
    Sc(i) = Viscosity/DensityMassGas/diffusivity0(i);
    Sh(i) = J_m*Re*Sc(i)^(1/3);
    k_mat(i) = Sh(i)*diffusivity0(i)/ParticleDiameter;
end

%% -- Governing Equations 
for i=1:NS
    pfr_het(i)=(3600*a_v_particle*k_mat(i)*DensityMassGas*(y(i+Nsolid)-y(i)))/G;
    pfr_het(i+Nsolid)=-(3600*a_v_particle*k_mat(i)*DensityMassGas*(y(i+Nsolid)-y(i)))+ ...
        NetRateProduction(i)*(1-VoidFraction)*CatalystDensity*mw(i);
end

RateH=0;
for j=1:NR
    RateH=RateH+deltaH(j)*ReactionRate(j); %kJ/kg_cat/h
end

pfr_het(Nsolid)=(3600*1e-3*h_thermal*a_v_particle*(TemperatureSurface-Temperature)+ ...
    U*(4/TubeDiameter)*(MoltenSaltsTemperature-Temperature))/G/SpecificHeatP;

pfr_het(NEQ-1)=(-RateH*(1-VoidFraction)*CatalystDensity - ...
    3600*1e-3*h_thermal*a_v_particle*(TemperatureSurface-Temperature))/G/SpecificHeatP;

pfr_het(NEQ)=-(((150*(1-VoidFraction)^2/VoidFraction^3))*...
    Viscosity*SuperficialVelocity/ParticleDiameter^2+...
    7/4*(DensityMassGas*SuperficialVelocity^2/ParticleDiameter)...
    *((1-VoidFraction)/VoidFraction^3))/1e5;

pfr_het=pfr_het';
end

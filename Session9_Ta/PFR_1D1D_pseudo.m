function PFR_1D1D_pseudo=PFR_1D1D_pseudo(t,y)

global Stoichiometry deltaH TubeDiameter ParticleDiameter ReactorLength...
    CatalystDensity VoidFraction MoltenSaltsTemperature...
    G SpecificHeatP Viscosity0 mw
global ExternalHeatCoefficient TubeThickness lambdaTube lambda0 lambda_cat
global diffusivity0 a_v_particle
global NS NEQ NR r_vector startSolidIndex epsi tau pore_diameter
global initial_point ParticleRadius

%molar_fraction:
mol=zeros(NS,1);
molSolid=zeros(NR,NS);
omega=zeros(NS,1);
omegaSolid=zeros(NR,NS);
MolarFractionSolid=zeros(NR,NS);
TemperatureSolid=zeros(NR,1);

for i=1:NS
    mol(i)=y(i)/mw(i); %kmol (assumed: bc=>1 kg)
    omega(i)=y(i);
     
    for j=1:NR
        molSolid(j,i)=y(startSolidIndex+(j-1)*NEQ+i)/mw(i);
        omegaSolid(j,i)=y(startSolidIndex+(j-1)*NEQ+i);
    end
end

MolarFraction=mol./(sum(mol));

for j=1:NR
    MolarFractionSolid(j,:)=molSolid(j,:)./(sum(molSolid(j,:)));
end

Temperature=y(NS+1);
for j=1:NR
    i=NEQ;
    TemperatureSolid(j)=y(startSolidIndex+(j-1)*NEQ+i);
end

Pressure_bar=y(end); %bar
Pressure_Pa=Pressure_bar*1e5; %Pa

DensityMolarGas=1e-3*Pressure_Pa/8.314/Temperature; %kmol/m3

mw_av=0;
for i=1:NS
    mw_av=mw_av+MolarFraction(i)*mw(i);
end

mw_solid=zeros(1,NR);
for j=1:NR
    for i=1:NS
        mw_solid(j)=mw_solid(j)+MolarFractionSolid(j,i)*mw(i);
    end
end

DensityMassGas=DensityMolarGas*mw_av; %kg/m3

SuperficialVelocity=G/DensityMassGas/3600; %m/s
%--
%KINETICS:
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
        Pressure_bar*MolarFractionSolid(j,2);

    for k=1:3                    
        for i=1:NS
            NetRateProduction(j,i)=NetRateProduction(j,i)+...
                Stoichiometry(k,i)*ReactionRate(k);
        end
        RateH(j)=RateH(j)+deltaH(k)*ReactionRate(k); %kJ/kg_cat/h
    end
end

eta = ones(NS,1);
r = (r_vector-initial_point)/(ParticleRadius-initial_point);
V=4/3*pi*r.^3;

for i = 1:NS
    Ri(i)=trapz(V,NetRateProduction(:,i)*CatalystDensity*mw(i)/3600)/(4/3*pi*(r(end)^3-r(1)^3));
    Rsurface(i)= NetRateProduction(end,i)*CatalystDensity*mw(i)/3600;
    eta(i)= Ri(i)/(Rsurface(i)+1e-12);
end

%--------
%% HEAT TRANSFER
Viscosity=Viscosity0;
lambda=lambda0;
diffusivity=diffusivity0;

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
U=U*3600/1000;

%% GOVERNING EQUATION
PFR_1D1D_pseudo=zeros(length(y),1);

for i=1:NS+1
    %Gas Phase
    if i<NS+1
        %Species Mass Balance
        PFR_1D1D_pseudo(i) = NetRateProduction(NR,i)*eta(i)*(1-VoidFraction)*CatalystDensity*mw(i)/G;
    else
        %Energy Balance with computed overall heat transfer coefficient
        PFR_1D1D_pseudo(NS+1) = (-RateH(NR)*(1-VoidFraction)*CatalystDensity + ...
    U*(4/TubeDiameter)...
    *(MoltenSaltsTemperature-Temperature))/G/SpecificHeatP;
    end
    
    D_eff_c=zeros(1,NR);
    D_eff_Kn=zeros(1,NR);
    D_eff=zeros(1,NR);
    if i<NS+1
        for j=1:NR
            %Bosanquet
            D_eff_c(j)=diffusivity(i)*epsi/tau;
            D_eff_Kn(j)=pore_diameter/3*(8*8.314*TemperatureSolid(j)/(pi*mw(i)*1e-3))^0.5;
            D_eff(j)=1/(1/D_eff_c(j)+1/D_eff_Kn(j));
        end
    end
    
    %Solid Species Balances
    for j=1:NR
        if j==1 %Catalyst Center
            if i<NS+1
                yp_r=(omegaSolid(j+1,i)*mw_solid(j+1)-omegaSolid(j,i)*mw_solid(j)) ...
                    /(r_vector(j+1)-r_vector(j));
                PFR_1D1D_pseudo(startSolidIndex+(j-1)*NEQ+i)=yp_r;
            else
                Tp_r=(TemperatureSolid(j+1)-TemperatureSolid(j))/(r_vector(j+1)-r_vector(j));
                PFR_1D1D_pseudo(startSolidIndex+(j-1)*NEQ+i)=Tp_r;
            end
        elseif(j>=2 && j<NR)
            if i<NS+1
                DoTp_r=(D_eff(j+1)/TemperatureSolid(j+1)-D_eff(j)/TemperatureSolid(j)) ...
                    /(r_vector(j+1)-r_vector(j));
                yp_r=(omegaSolid(j+1,i)*mw_solid(j+1)-omegaSolid(j,i)*mw_solid(j)) ...
                    /(r_vector(j+1)-r_vector(j));
                ypp_r=(omegaSolid(j+1,i)*mw_solid(j+1)-2*omegaSolid(j,i)*mw_solid(j)+omegaSolid(j-1,i)*mw_solid(j-1)) ...
                    /((r_vector(j+1)-r_vector(j))*(r_vector(j)-r_vector(j-1))); 
                diffusionTerm = (1/r_vector(j)*Pressure_Pa*1e-3/8.3144) ...
                    *(2*D_eff(j)/TemperatureSolid(j)*yp_r+r_vector(j)*(DoTp_r*yp_r+D_eff(j)/TemperatureSolid(j)*ypp_r));
                reactionTerm = NetRateProduction(j,i)*CatalystDensity*mw(i)/3600;
                PFR_1D1D_pseudo(startSolidIndex+(j-1)*NEQ+i)=(diffusionTerm+reactionTerm);
            else
                Tp_r=(TemperatureSolid(j+1)-TemperatureSolid(j))/(r_vector(j+1)-r_vector(j));
                Tpp_r=(TemperatureSolid(j+1)-2*TemperatureSolid(j)+TemperatureSolid(j-1))/((r_vector(j+1)-r_vector(j))*(r_vector(j)-r_vector(j-1))); 
                PFR_1D1D_pseudo(startSolidIndex+(j-1)*NEQ+i)=(lambda_cat/r_vector(j)*(r_vector(j)*Tpp_r+2*Tp_r)-RateH(j)*CatalystDensity/3600*1000);
            end
        elseif j==NR %Gas-Solid Interface
            if i<NS+1           
                 %Species Balance Equation
                  PFR_1D1D_pseudo(startSolidIndex+(j-1)*NEQ+i)= y(startSolidIndex+(j-1)*NEQ+i)-y(i);
            else
                 %Energy Balance Equation
                  PFR_1D1D_pseudo(startSolidIndex+(j-1)*NEQ+i)= y(startSolidIndex+(j-1)*NEQ+i)-y(i);
            end
        end
    end
end

%--Ergun Equation:
PFR_1D1D_pseudo(end)=-(150*(1-VoidFraction)^2/(VoidFraction^3)...
    *Viscosity*SuperficialVelocity/ParticleDiameter^2+...
    +7/4*DensityMassGas*SuperficialVelocity^2/ParticleDiameter*...
    (1-VoidFraction)/(VoidFraction^3))/1e5;

end



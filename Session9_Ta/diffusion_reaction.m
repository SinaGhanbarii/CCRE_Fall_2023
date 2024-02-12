function diffusion_reaction=diffusion_reaction(y)

%
% CCRE - Prof. Matteo Maestri - November 20 - a.a 2023/24
%

global NR NS NEQ mw N_reactions deltaH Stoichiometry diffusivity0 epsi tau pore_diameter ...
    delta_r CatalystDensity pressure r_vector lambda_cat yIN

% i => index for the species
% j => index for the points inside the domain

molSolid=zeros(NR,NS);
omegaSolid=zeros(NR,NS);
MolarFractionSolid=zeros(NR,NS);
TemperatureSolid=zeros(NR,1);

for i=1:NS
    for j=1:NR
        molSolid(j,i)=y((j-1)*NEQ+i)/mw(i);
        omegaSolid(j,i)=y((j-1)*NEQ+i);
    end
end

for j=1:NR
    MolarFractionSolid(j,:)=molSolid(j,:)./(sum(molSolid(j,:)));
end

for j=1:NR
    i=NEQ;
    TemperatureSolid(j)=y((j-1)*NEQ+i);
end

mw_solid=zeros(1,NR);
for j=1:NR
    for i=1:NS
        mw_solid(j)=mw_solid(j)+MolarFractionSolid(j,i)*mw(i);
    end
end

%-KINETIC SCHEME:
NetRateProduction=zeros(NR,NS);
RateH=zeros(NR,1);
Pressure_bar = pressure/1e5; % bar

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

%% -- Governing Equations
D_eff_c=zeros(NR,NS);
D_eff_Kn=zeros(NR,NS);
D_eff=zeros(NR,NS);

for i=1:NEQ
    if i<NEQ
        for j=1:NR
            %Bosanquet
            D_eff_c(j,i)=diffusivity0(i)*epsi/tau;
            D_eff_Kn(j,i)=pore_diameter/3*(8*8.314*TemperatureSolid(j)/(pi*mw(i)*1e-3))^0.5;
            D_eff(j,i)=1/(1/D_eff_c(j,i)+1/D_eff_Kn(j,i));
        end
    end
    
    %Solid Species Balances
    for j=1:NR
        if j==1 %Catalyst Center
            if i<NEQ
                % species
                yp_r=(omegaSolid(j+1,i)*mw_solid(j+1)-omegaSolid(j,i)*mw_solid(j))/delta_r;
                diffusion_reaction((j-1)*NEQ+i)=yp_r;
            else
                % energy
                Tp_r=(TemperatureSolid(j+1)-TemperatureSolid(j))/delta_r;
                diffusion_reaction((j-1)*NEQ+i)=Tp_r;
            end
        elseif(j>=2 && j<NR)
            if i<NEQ
                % species
                DoTp_r=(D_eff(j+1,i)/TemperatureSolid(j+1)-D_eff(j,i)/TemperatureSolid(j))/delta_r;
                yp_r=(omegaSolid(j+1,i)*mw_solid(j+1)-omegaSolid(j,i)*mw_solid(j))/delta_r;
                ypp_r=(omegaSolid(j+1,i)*mw_solid(j+1)-2*omegaSolid(j,i)*mw_solid(j)+ ...
                    omegaSolid(j-1,i)*mw_solid(j-1))/delta_r^2; 
                diffusionTerm = 1/r_vector(j)^2*pressure*1e-3/8.3144*(2*r_vector(j)*D_eff(j,i) ...
                    /TemperatureSolid(j)*yp_r+r_vector(j)^2*(DoTp_r*yp_r+D_eff(j,i) ...
                    /TemperatureSolid(j)*ypp_r));
                reactionTerm = NetRateProduction(j,i)*CatalystDensity*mw(i)/3600;
                diffusion_reaction((j-1)*NEQ+i)=(diffusionTerm+reactionTerm);
            else
                % energy
                Tp_r=(TemperatureSolid(j+1)-TemperatureSolid(j))/delta_r;
                Tpp_r=(TemperatureSolid(j+1)-2*TemperatureSolid(j)+TemperatureSolid(j-1))/delta_r^2; 
                diffusion_reaction((j-1)*NEQ+i)=(lambda_cat/r_vector(j)*(2*Tp_r+r_vector(j)*Tpp_r) ...
                    -RateH(j)*CatalystDensity/3600*1000);
            end
        else %Gas-Solid Interface
            if i<NEQ             
                diffusion_reaction((j-1)*NEQ+i)=y((j-1)*NEQ+i)-yIN(i);  
            else
                diffusion_reaction((j-1)*NEQ+i)=y((j-1)*NEQ+i)-yIN(end);
            end
        end
    end
end

diffusion_reaction=diffusion_reaction';
end


function fluidized_bed=fluidized_model(x)


% CCRE - Prof. Matteo Maestri - November 27 - a.a 2023/24
%
% KUNII-LEVENSPIEL MODEL
% BUBBLING BED REACTOR (GELDART A)
%

global NS NP L NTOT delta_z z_star v_bubble v_mf omega_0 v_0 delta integral ...
    yp_bubble rho_0 MW epsi_mf v_single_bubble d_bubbles g diffusivity ...
    k_bc k_ce P0 T0 catalyst_density Rate_emulsion Rate_bubble Rate_cloud ...
    gamma_b gamma_c gamma_e f_cl f_w

%% RIARRANGE DATA IN MATRIX FORM
y_bubble=zeros(NP,NS);
yp_bubble=zeros(NP,NS);
y_cloud=zeros(NP,NS);
y_emulsion=zeros(NP,NS);
integral_function=zeros(NP,NS);
integral=zeros(NS,1);
k_bc=zeros(NP,NS);
k_ce=zeros(NP,NS);
rho_bubble=zeros(NP,1);
rho_cloud=zeros(NP,1);
rho_emulsion=zeros(NP,1);

fluidized_bed=zeros(size(x));

for i=1:NP
    for k=1:NS
        j=1;
            y_bubble(i,k)=x(NTOT*(j-1)+NS*(i-1)+k);
        j=2;
            y_cloud(i,k)=x(NTOT*(j-1)+NS*(i-1)+k);
        j=3;
            y_emulsion(i,k)=x(NTOT*(j-1)+NS*(i-1)+k);
    end
end

%% DERIVATIVES COMPUTATION (ONLY FOR BUBBLE PHASE)
for k=1:NS
    %inner points:
    for i=2:NP-1
        yp_bubble(i,k)=(y_bubble(i,k)-y_bubble(i-1,k))/(z_star(i)-z_star(i-1));
    end

    %initial point (reactor inlet):
    i=1;
        yp_bubble(i,k)=(y_bubble(i+1,k)-y_bubble(i,k))/(z_star(i+1)-z_star(i));

    %final point (fluidized bed outlet):
    i=NP;
        yp_bubble(i,k)=(y_bubble(i,k)-y_bubble(i-1,k))/(z_star(i)-z_star(i-1));
end

%% PROPERTIES COMPUTATION
for i=1:NP
    rho_bubble(i)=rho_0;
    rho_cloud(i)=rho_0;
    rho_emulsion(i)=rho_0;
end

for i=1:NP
    for k=1:NS
        k_bc(i,k) = 4.*(v_mf/d_bubbles)+5.85*(diffusivity^0.5*g^0.25/d_bubbles^(5/4));
        k_ce(i,k) = 6.77*(diffusivity*epsi_mf*v_bubble/d_bubbles^3)^0.5;
    end
end

%% COMPUTATION OF THE INTEGRAL FOR THE EMULSION BALANCE
for i=1:NP
    for k=1:NS
        rho_interface=(rho_cloud(i)+rho_emulsion(i))/2;
        integral_function(i,k)= k_ce(i,k)*rho_interface*(y_cloud(i,k)-y_emulsion(i,k));
    end
end

for k=1:NS
    integral(k)= trapz(z_star, integral_function(:,k));
end

%% COMPUTATION OF KINETICS
%n-butane partial oxidation to maleic anhydride
%reaction 1: n-butane + 3.5O2 -> maleic anhydride(MA) + 4H2O (MA synthesis)
%reaction 2: n-butane + 5.5O2 -> 2CO +2CO2 + 5H2O (n-butane oxidation)
%reaction 3: MA + O2 -> 4CO + H2O (MA partial oxidation)

%species order: 1)n-butane 2)O2 3)MA 4)CO 5)CO2 6)H2O 7)N2
nu=[-2 1 1 0];

Rate_bubble=zeros(NP,NS);
Rate_cloud=zeros(NP,NS);
Rate_emulsion=zeros(NS,1);
xi_bubble=zeros(NP,NS);
xi_cloud=zeros(NP,NS);
xi_emulsion=zeros(NP,NS);

for i=1:NP
    
    %conversion of mass fraction in molar fraction
    sum_bubble=0;
    sum_cloud=0;
    sum_emulsion=0;
    for k=1:NS
       molarfraction_bubble(k)=y_bubble(i,k)/MW(k);
       molarfraction_cloud(k)=y_cloud(i,k)/MW(k);
       molarfraction_emulsion(k)=y_emulsion(i,k)/MW(k);
       sum_bubble=sum_bubble+molarfraction_bubble(k);
       sum_cloud=sum_cloud+molarfraction_cloud(k);
       sum_emulsion=sum_emulsion+molarfraction_emulsion(k);
    end
    molarfraction_bubble=molarfraction_bubble/sum_bubble;
    molarfraction_cloud=molarfraction_cloud/sum_cloud;
    molarfraction_emulsion=molarfraction_emulsion/sum_emulsion;
    
    for k=1:NS
        xi_bubble(i,k)=max(molarfraction_bubble(k),0);
        xi_cloud(i,k)=max(molarfraction_cloud(k),0);
        xi_emulsion(i,k)=max(molarfraction_emulsion(k),0);
    end
       
    %computation of reaction rates parameters

    Trif=273.15+288; %[K]
    k_r = 21900;
    K_eq = 11.5;
    K_MeOH = 14000;
    K_H2O = 47000;
    Rgas=8.314/1e3; %[kJ/mol]
    
    disp(xi_bubble(i,1))
    %computation of reaction rates
    C_tot = P0/Rgas/Trif;
    disp(C_tot)
    r_Bubble = (k_r*(C_tot^2*xi_bubble(i,1)^2))*(1-(C_tot^2*xi_bubble(i,2)*xi_bubble(i,3))/(C_tot^2*xi_bubble(i,1)^2*K_eq))/(1+K_MeOH*xi_bubble(i,1)*C_tot+K_H2O*C_tot*xi_bubble(i,3)).^2;
    r_cloud = (k_r*(C_tot^2*xi_cloud(i,1)^2))*(1-(C_tot^2*xi_cloud(i,2)*xi_cloud(i,3))/(C_tot^2*xi_cloud(i,1)^2*K_eq))/(1+K_MeOH*xi_cloud(i,1)*C_tot+K_H2O*C_tot*xi_cloud(i,3)).^2;
    r_emulsion = (k_r*(C_tot^2*xi_emulsion(i,1)^2))*(1-(C_tot^2*xi_emulsion(i,2)*xi_emulsion(i,3))/(C_tot^2*xi_emulsion(i,1)^2*K_eq))/(1+K_MeOH*xi_emulsion(i,1)*C_tot+K_H2O*C_tot*xi_emulsion(i,3)).^2;
    
    %computation and storing of species molar production rates
    rk_prodution_bubble=nu*r_Bubble; %kmol/kgcat/s
    rk_prodution_cloud=nu*r_cloud; %kmol/kgcat/s
    rk_prodution_emulsion=nu*r_emulsion; %kmol/kgcat/s
    
    for k=1:NS
        Rate_bubble(i,k)=rk_prodution_bubble(k)*catalyst_density; %kmol/m3cat/s
        Rate_cloud(i,k)=rk_prodution_cloud(k)*catalyst_density; %kmol/m3cat/s
        Rate_emulsion(i,k)=rk_prodution_emulsion(k)*catalyst_density; %kmol/m3cat/s
    end
end

%% GOVERNING EQUATIONS 
for i=1:NP
    rho_interface_bc=(rho_bubble(i)+rho_cloud(i))/2;
    rho_interface_ce=(rho_cloud(i)+rho_emulsion(i))/2;
    for k=1:NS
        j=1; % bubble phase
            if i == 1
                fluidized_bed(NTOT*(j-1)+NS*(i-1)+k)=y_bubble(i,k)-omega_0(k);
            else
                fluidized_bed(NTOT*(j-1)+NS*(i-1)+k)= -rho_0*v_bubble*yp_bubble(i,k) ...
                +k_bc(i,k)*rho_interface_bc*(y_cloud(i,k)-y_bubble(i,k)) ...
                +Rate_bubble(i,k)*MW(k)*gamma_b;
            end
        j=2; % cloud phase
            fluidized_bed(NTOT*(j-1)+NS*(i-1)+k)= -k_bc(i,k)*rho_interface_bc*(y_cloud(i,k)-y_bubble(i,k)) ...
                +k_ce(i,k)*rho_interface_ce*(y_emulsion(i,k)-y_cloud(i,k)) ...
                +Rate_cloud(i,k)*MW(k)*gamma_c;
       
        j=3; % emulsion phase
            fluidized_bed(NTOT*(j-1)+NS*(i-1)+k)= Rate_emulsion(i,k)*MW(k)*gamma_e ...
                +k_ce(i,k)*rho_interface_ce*(y_cloud(i,k)-y_emulsion(i,k)); 
    end
end
end


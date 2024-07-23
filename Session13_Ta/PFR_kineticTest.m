% 
% - CCRE - PROF. MATTEO MAESTRI - 2023-2024 - December 11
%

clear all, clc 

global kr epsi

k0 = 1.96e9;
Ea_R = 8.491e3;
epsi = 0.4;

Pressure=101325; %Pa
R=8.314; %J/mol/K
Q_298K_1atm=2e-2; %Nm3/s
d_r=0.02; %m reactor diameter
L_r=0.06; %m reactor length
A_r=pi*d_r^2/4; %m2
V_r=A_r*L_r; %m3

c_0= [0.75 0.1]; %kmol/m3

Temperature=550; %K
Q=Q_298K_1atm*(Temperature/298); %m3/s
tau=V_r/Q; %s

kr=k0*exp(-Ea_R/Temperature);

[t,y]=ode15s(@PFR,[0 tau],c_0);

conversion_A=(y(1,1)-y(end,1))/y(1,1)*100;
CAout=y(end,1);
CBout=y(end,2);

disp("temperatura: "+Temperature)
disp("CAout: "+CAout)
disp("CBout: "+CBout)
disp("conversion_A: "+conversion_A+" %")

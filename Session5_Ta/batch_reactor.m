function batch_reactor = batch_reactor(t,y)
global E1 R k1_0 E2 k2_0 deltaH1 deltaH2 U A_v rho cp t_s

cA = y(1);
cP = y(2);
cC = y(3);
T = y(4);

k1 = k1_0*exp(-E1/R/T);
k2 = k2_0*exp(-E2/R/T);
r1 = k1*cA; % mol/m3/s
r2= k2*cP;

if t<t_s
    T_h = 345; %K
else 
    T_h = 295; %K
end

batch_reactor(1)=-r1;
batch_reactor(2)=r1-r2;
batch_reactor(3)=r2;
batch_reactor(4)=(-r1*deltaH1-r2*deltaH2+U*A_v*(T_h - T))/rho/cp;

batch_reactor = batch_reactor';

end
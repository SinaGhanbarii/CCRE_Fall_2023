function semi_batch_reactor = semi_batch_reactor(t,y)
global E1 R k1_0 E2 k2_0 deltaH1 deltaH2 U A_v rho cp t_s Qin TFeed cAin t_d

cA = y(1);
cP = y(2);
cC = y(3);
T = y(4);
V = y(5);

k1 = k1_0*exp(-E1/R/T);
k2 = k2_0*exp(-E2/R/T);
r1 = k1*cA; % mol/m3/s
r2= k2*cP;

if t<t_s
    T_h = 345; %K
else 
    T_h = 295; %K
end

if t < t_d
    semi_batch_reactor(1)=-r1;
    semi_batch_reactor(2)=r1-r2;
    semi_batch_reactor(3)=r2;
    semi_batch_reactor(4)=(-r1*deltaH1-r2*deltaH2+U*A_v*(T_h-T))/rho/cp;
    semi_batch_reactor(5)=0;
else
    semi_batch_reactor(1)=-r1 + (cAin-cA)*Qin/V;
    semi_batch_reactor(2)=r1-r2-cP*Qin/V;
    semi_batch_reactor(3)=r2-cC*Qin/V;
    semi_batch_reactor(4)=(-r1*deltaH1-r2*deltaH2+U*A_v*(T_h-T)+rho*cp*Qin*(TFeed-T)/V)/rho/cp;
    semi_batch_reactor(5)=Qin;
end

semi_batch_reactor = semi_batch_reactor';

end
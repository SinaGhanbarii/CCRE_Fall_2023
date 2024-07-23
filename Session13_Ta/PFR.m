function PFR=PFR(t,y)
global kr epsi

cA = y(1);
rate = kr*cA*(1-epsi);

PFR(1)=-rate;
PFR(2)=rate;

PFR=PFR';
end
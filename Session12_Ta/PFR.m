function PFR=PFR(t,y)
global kr

cA = y(1);
rate = kr*cA;

PFR(1)=-rate;
PFR(2)=rate;

PFR=PFR';
end
function exe1_adi_tubular = exe1_adi_tubular(V,y)
global NS nu index DHr_ref Cpmix R Pressure

T = y(NS+1);
Tref = 25+273;

ntot_dot = 0;
for i = 1:NS
    ntot_dot = ntot_dot + y(i);
end

for i = 1:NS
    x(i) = y(i)/ntot_dot;
end

ctot = Pressure*1e-3/R/T;

for i = 1:NS
    conc(i) = x(i)*ctot;
end

DHr = DHr_ref + Cpmix*(T-Tref);

% Coefficients for Kinetic Constant dependency on temperature!
a = 0.3099;
b = -275.06;
c = 68248;
% Assign the value of k_r using the given Temp.
% kr = polyval([a b c], T);
kr = a*T^2+b*T+c;
rate = kr * conc(index('O2')) * conc(index('NO'))^2;
for i = 1:NS
    exe1_adi_tubular(i) = (nu(i) * rate);
end

i = NS+1;
exe1_adi_tubular(i) = -DHr*rate/(Cpmix*ntot_dot);

exe1_adi_tubular = exe1_adi_tubular';
end
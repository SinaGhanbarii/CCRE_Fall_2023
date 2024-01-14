function exe1_iso_tubular = exe1_iso_tubular(V,n)
global NS nu kr index ctot

ntot_dot = sum(n);
x = n/ntot_dot;
c = x*ctot;

kr = 1.4e4;
rate = kr* c(index('O2')) * c(index('NO'))^2;

for i= 1:NS
    exe1_iso_tubular(i) = (nu(i)*rate);
end
exe1_iso_tubular = exe1_iso_tubular';

end
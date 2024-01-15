function pfr_2=pfr_2(t,y)
global k_r nu 
global alpha1 alpha2 alpha3

%this is to ensure conc of A should not become negative 
y(1) = max(y(1),1e-6);

%reaction rates 
r(1)=k_r(1)*y(1)^alpha1;
r(2)=k_r(2)*y(1)^alpha2;
r(3)=k_r(3)*y(1)^alpha3;

R=zeros(4,1);
for i=1:4
    R(i)=0;
    for j=1:3
        R(i)=R(i)+nu(j,i)*r(j);
    end
end
for i=1:4
    pfr_2(i)=R(i);
end
pfr_2=pfr_2';
end
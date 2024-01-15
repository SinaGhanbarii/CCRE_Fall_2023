function pfr = pfr(t,y)
global k_r nu
global alpha1 alpha2

% This is to ensure that the A isn't negative
y(1) = max(y(1), 1e-6);

%Reaction rates
r(1) = k_r(1)*y(1)^alpha1;
r(2) = k_r(2)*y(1)^alpha2;

R = ones(3,1);

for i = 1:3
    R(i) = 0;
    for j=1:2
        R(i) = R(i) + nu(j,i)*r(j);
    end
end

for i=1:3
    pfr(i) = R(i);
end

pfr = pfr';
end
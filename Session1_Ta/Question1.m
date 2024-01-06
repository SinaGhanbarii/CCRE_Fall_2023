%Excercise 1
clear, clc
global sigma_tetta_square
time = [0,0.15,0.29,0.44,0.59,0.74,0.88,1.03,1.18,1.32,1.47,1.62,1.76,1.91, ...
    2.06,2.21,2.35,2.50]; % [s]

C = [0.8364,1.3283,1.9911,2.8169,3.7614,4.7402,5.6382,6.3295,6.7063,6.7063,6.3295,	...
    5.6382,4.7402,3.7614,2.8170,1.9911,1.3283,0.8364]; % [kmol/m3]

figure
plot(time,C, 'o')
xlabel('Time (s)')
ylabel('Concentration (kmol/m3)')

%%Data
d = 0.04; %m
A = pi*0.25*d^2; %m2
L = 30; %m
v = 8.1; %m/s
V = A*L; %m3
Q = A*v; %m3/s

kr = 0.1; %m3/kmol/s

%% Getting dispersion Coeff.

C_integral = trapz(time,C);
N = length(C);

%C_integral2 = 0;

%for i=1:N-1
%    C_integral2 = C_integral2 + (C(i)+C(i+1))/2 * (time(i+1)-time(i));
%end

E = C/C_integral;

%zero order moment

m0 = trapz(time,E);

%RTD
W = E/m0;

%1st order moment
E_time = ones(1,length(E));
for i= 1:length(E)
    E_time(i) = E(i)*time(i);
end

t_mean = trapz(time,E_time);

%2nd order moment
E_time_square = ones(1, length(E));
for i=1:length(E)
    E_time_square(i) = E(i)* time(i)^2;
end

sigma_square_integral = trapz(time, E_time_square);
sigma_square = sigma_square_integral - t_mean^2;
sigma_tetta_square = sigma_square/t_mean^2 ;

first_guess = 2/sigma_tetta_square;

Pe = fsolve(@PecletFunction, first_guess);
D_eff = v*L/Pe; 
fprintf('D_eff : %f m2/s', D_eff)

%% Function
function PecletFunction = PecletFunction(x)
global sigma_tetta_square
Pe = x;

PecletFunction = Pe - 2/sigma_tetta_square*(1-1/Pe*(1-exp(-Pe)));
end
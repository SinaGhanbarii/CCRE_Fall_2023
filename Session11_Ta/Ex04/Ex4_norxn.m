close all; clear all; clc;

global b c DiffA DiffB HA KLa KGa

%% Data struct
b = 1;
c = 1;
DiffA = 1e-6; % diffusion coefficient [m2/hr]
DiffB = DiffA; % diffusion coefficient [m2/hr]

HA = 12.5; % [Pa m3/mol]

KLa = 0.1; % mass transfer in the liquid phase [1/hr]
KGa = 0.32; % mass transfer in the gas phase [mol/m3/hr/Pa]

%% Input data
pTot = 100000; % [Pa]
Ctot = 56000; % [mol/m3]
pAin = 100; % [Pa]
CAin = 0; %[mol/m3]

FtotG = 10^5; %[mol/hr/m2]
FtotL = 7*10^5; %[mol/hr/m2]

for j = 500:0.5:520
    Lreactor = j;
    %% Shooting algorithm
    maxinterations = 200; % iterations
    alpha = 0.1;
    deltaPthreshold = 0.1; % [Pa]
        
    pAout = (1-0.8)*pAin; % first guess solution (from co-current)
        
    for i = 1:maxinterations
           %ODE solution
           Yin = [pAout,CAin]';
           [z , Y] = ode45(@TowerCounterCurrent,[0 Lreactor], Yin, [],pTot,Ctot,FtotG,FtotL);
        
           pA = Y(:,1);
        
           %check solution
           error = pA(end) - pAin;
        
           if (abs(error) < deltaPthreshold)
                break;
           end
        
           %new guess
           pAout = pAout - alpha*error;
    end
    pA = Y(:,1);    % partial pressure of A [Pa]
    cA = Y(:,2);    % concentration of B [kmol/m3]
    fprintf("Reactor length is : %f  [m] \n", Lreactor)
    fprintf("outlet pressure is : %f  [Pa] \n", pA(1))
    if pA(1) < 20
        break
    end 
end

%% Post-processing
figure;
plot(z,pA);
hold on;
plot(z,cA); 
xlabel('axial coordinate [m]');
ylabel('partial pressure of A [Pa]');


%% ODE system (defined locally, requires at least MATLAB-2016b)
function dY = TowerCounterCurrent(z,Y, pTot, Ctot, FtotG, FtotL)
global b c DiffA DiffB HA KLa KGa

pA = Y(1);
CA = Y(2);

r  = OverallRateOfChange_norxn(pA,CA);

dpA_over_dz = pTot/FtotG*r;
dCA_over_dz = Ctot/FtotL*r;

dY = [dpA_over_dz,dCA_over_dz]'; 

end
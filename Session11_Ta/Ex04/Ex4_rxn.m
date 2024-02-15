close all; clear variables;

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
CBin = 800; %[mol/m3]

FtotG = 10^5; %[mol/hr/m2]
FtotL = 7*10^5; %[mol/hr/m2]

for j = 1:0.1:10
    Lreactor = j;
    %% Shooting algorithm
    maxinterations = 100; % iterations
    alpha = 0.1;
    deltaPthreshold = 0.1; % [Pa]
    
    pAout = (1-0.8)*pAin; % first guess solution (from co-current)
    
    for i = 1:maxinterations
    
           % ODE solution
           Yin = [pAout,CBin]';
           [z , Y] = ode45(@TowerCounterCurrent,[0 Lreactor], Yin, [],pTot,Ctot,FtotG,FtotL);
    
           pA = Y(:,1);
    
           % check solution
           error = pA(end) - pAin;
    
           if (abs(error) < deltaPthreshold)
                break;
           end
    
           % new guess
           pAout = pAout - alpha*error;
    end
    
    pA = Y(:,1);    % partial pressure of A [Pa]
    CB = Y(:,2);    % concentration of B [kmol/m3]
    if pA(1) < 20
        fprintf("Reactor lengtrh is : %f [m] \n",  Lreactor)
        break
    end 
end

%% Post-processing
figure;
xlabel('axial coordinate [m]');
hold on;
ylabel('partial pressure of A [Pa]'); 
plot(z,pA); 

figure;
xlabel('axial coordinate [m]');
hold on;
yyaxis left; ylabel('partial pressure of A [Pa]'); 
plot(z,pA); 
yyaxis right; ylabel('concentration of B [mol/m3]');
plot(z,CB);

%% Reconstructing the overall rate of change and the enhancing factor
for i=1:size(pA,1)
    [r(i),E(i), resG(i), resL(i)] = OverallRateOfChange_rxn(pA(i), CB(i));
end

figure;
plot(z,resG);
hold on 
plot(z,resL);
xlabel('axial coordinate [m]');
ylabel('series resistances'); 
legend('Gas','Film')

for i=1:size(pA,1)
    pAinterface(i) = pA(i) - r(i)*resG(i);
    cAinterface(i) = pAinterface(i)/HA;
    cAbulk(i) = cAinterface(i) - r(i)*(resL(i)/HA);
end 

figure
plot(z,pA);
hold on 
plot(z,pAinterface);
xlabel('axial coordinate [m]');
ylabel('series conc [mol/m3]'); 
legend('Gas','Interface')

%% ODE system (defined locally, requires at least MATLAB-2016b)
function dY = TowerCounterCurrent(z,Y, pTot, Ctot, FtotG, FtotL)
global b c DiffA DiffB HA KLa KGa

pA = Y(1);
CB = Y(2);

r = OverallRateOfChange_rxn(pA,CB);

dpA_over_dz = pTot/FtotG*r;
dCB_over_dz = -Ctot/FtotL*b*r;

dY = [dpA_over_dz,dCB_over_dz]'; 

end
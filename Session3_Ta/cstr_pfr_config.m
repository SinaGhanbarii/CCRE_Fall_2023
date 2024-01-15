clc
clear all
close all
hold off
global k_r nu
global Q V tau_total tau_cstr tau_pfr c_0
global alpha1 alpha2 alpha3

%reaction orders for parallel reactions 
alpha1 = 2; %B
alpha2 = 1; %C
alpha3 = 0.5; %D

%rate constants for parallel reactions 
k_r=[1 1 1];

%total reactor volume 
V=10; %m3

%stoichiometry matrix
nu=[-1 1 0 0
    -1 0 1 0
    -1 0 0 1];

%initial concentration of reactant in feed
c_0_in=[10 0 0 0];  %mol/m3

%initial concentration at inlet of a reactor
c_0=c_0_in;

%flow rate 
Q=10; %m3/h

%space time for the tubular flow reactor 
tau_total = V/Q;

i=0;

%%%%%%%%%%%%%%%%% CSTR + PFR configuration %%%%%%%%%%%%%%%%%%%%%%%%%%

for frac = 0.001:0.001:0.999
    i = i+1;
    fraction(i) = frac;

    tau_cstr = tau_total * frac;
    tau_pfr  = tau_total * (1- frac);
    
    %cstr + pfr :
    c_0 = c_0_in;
    y_cstr_2=fsolve('cstr_2',c_0);
    [t,y_pfr_2]=ode23s('pfr_2',[0 tau_pfr],y_cstr_2);
    
    conversion_A(i)=(c_0(1)-y_pfr_2(end,1))/c_0(1);
    
    selectivity_C(i)=y_pfr_2(end,3)/(c_0(1)-y_pfr_2(end,1));
    selectivity_B(i)=y_pfr_2(end,2)/(c_0(1)-y_pfr_2(end,1));
    selectivity_D(i)=y_pfr_2(end,4)/(c_0(1)-y_pfr_2(end,1));

    yield_C_cstr_pfr(i)=conversion_A(i)*selectivity_C(i);
    yield_B_cstr_pfr(i)=conversion_A(i)*selectivity_B(i);
    yield_D_cstr_pfr(i)=conversion_A(i)*selectivity_D(i); 

end

line_color = ['k' 'b' 'g' 'r'];
figure(1);
grid on;
title('CSTR + PFR')
hold on 
%plot(fraction, conversion_A, "Color",line_color(1), 'LineWidth',3)
%plot(fraction, yield_B_cstr_pfr, "Color",line_color(2), 'LineWidth',3)
plot(fraction, yield_C_cstr_pfr, "Color",line_color(3), 'LineWidth',3)
%plot(fraction, yield_D_cstr_pfr, "Color",line_color(4), 'LineWidth',3)
xlabel('frac (PFR ------> CSTR) ');
ylabel('C yield');
%legend('A','B', 'C', 'D')

%finding maximum conversion and the corresponding fraction of volume split 
[val, idx] = max(yield_C_cstr_pfr);
max_yield = val;
optimal_fraction = fraction(idx) ; 
fprintf('The optimal C yield is : %f \n', max_yield)
fprintf('The optimal fraction is : %f \n', optimal_fraction)

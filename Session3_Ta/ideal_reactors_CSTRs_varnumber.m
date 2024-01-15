clc
clear all
close all

global k_r nu
global Q V tau_total tau_cstr tau_pfr c_0 
global alpha1 alpha2

%reaction orders for parallel reactions 
alpha1 = 2; %B
alpha2 = 1; %C

%rate constants for parallel reactions 
k_r=[1 1];

%reactor volume  
V=10; %m3

%stoichiometry matrix
nu=[-1 1 0 
    -1 0 1 ];

%initial concentration of reactant in feed
c_0_in=[10 0 0];  %mol/m3

%flow rate 
Q=10; %m3/h
%space time for the tubular flow reactor 
tau_total_initial = V/Q;
tau_total = tau_total_initial;

%%%%%%%%%%%%%%%%% CSTR %%%%%%%%%%%%%%%%%%%%%%%%%%
%initial concentration at inlet of a reactor
c_0=c_0_in;

%cstr:
y_cstr=fsolve('cstr',c_0);
conversion_A_cstr=(c_0(1)-y_cstr(1))/c_0(1);
selectivity_B_cstr=y_cstr(2)/(c_0(1)-y_cstr(1));
selectivity_C_cstr=y_cstr(3)/(c_0(1)-y_cstr(1));
yield_B_cstr=conversion_A_cstr*selectivity_B_cstr;
yield_C_cstr=conversion_A_cstr*selectivity_C_cstr;

%%%%%%%%%%%%%%%%% PFR %%%%%%%%%%%%%%%%%%%%%%%%%%
%initial concentration at inlet of a reactor
c_0=c_0_in;

%pfr:
[t,y_pfr]=ode23s('pfr',[0 tau_total],c_0);
conversion_A_pfr=(y_pfr(1,1)-y_pfr(end,1))/y_pfr(1,1);
selectivity_B_pfr=y_pfr(end,2)/(y_pfr(1,1)-y_pfr(end,1));
selectivity_C_pfr=y_pfr(end,3)/(y_pfr(1,1)-y_pfr(end,1));
yield_B_pfr=conversion_A_pfr*selectivity_B_pfr;
yield_C_pfr=conversion_A_pfr*selectivity_C_pfr;

%%%%%%%%%%%%%%%%% CSTR in series %%%%%%%%%%%%%%%%%%%%%%%%%%
%no of CSTRs in series 
n_cstr_seq = [1, 3, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100]; 

%cstrs in series:

for count = 1:length(n_cstr_seq)
    %initial concentration at inlet of a reactor
    c_0=c_0_in;
    n_cstr = n_cstr_seq(count);
    tau_total = tau_total_initial/n_cstr;

    i=0;
    for no = 1:n_cstr
        i = i+1;
        y_cstr_frac=fsolve('cstr',c_0);
        y_cstr_seq(i,:) = y_cstr_frac;
        c_0=y_cstr_frac;
    end 

    conversion_A_loop=(c_0_in(1)-y_cstr_seq(end, 1))/c_0_in(1);
    selectivity_B_loop=y_cstr_seq(end,2)/(c_0_in(1)-y_cstr_seq(end,1));
    selectivity_C_loop=y_cstr_seq(end,3)/(c_0_in(1)-y_cstr_seq(end,1));
    yield_B_cstr_loop=conversion_A_loop*selectivity_B_loop;
    yield_C_cstr_loop=conversion_A_loop*selectivity_C_loop;

    selectivity_B_count(count) = selectivity_B_loop;
    yield_B_count(count) = yield_B_cstr_loop;
end


figure(1);
title('CSTR in series')
hold on
yline(selectivity_B_cstr,'-','CSTR')
yline(selectivity_B_pfr,'-','PFR')
plot(n_cstr_seq, selectivity_B_count, "Color",'k', 'LineWidth',3)
xlabel('increasing no of CSTRs');
ylabel('B selectivity');

figure(2);
title('CSTR in series')
hold on
yline(yield_B_cstr,'-','CSTR')
yline(yield_B_pfr,'-','PFR')
plot(n_cstr_seq, yield_B_count, "Color",'k', 'LineWidth',3)
xlabel('increasing no of CSTRs');
ylabel('B yield');



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


% Compute the space time

tau_total = V/Q;

%% CSTR
% Initial Concentration at inlet of reactor
c_0 = c_0_in;
%cstr 
y_cstr = fsolve('cstr', c_0);
conversion_A_cstr = (c_0(1) - y_cstr(1))/c_0(1);
Selectivity_B_cstr = y_cstr(2)/(c_0(1) - y_cstr(1));
yield_B_cstr = conversion_A_cstr * Selectivity_B_cstr;

%% PFR
% Initial Concentration at inlet of reactor

c_0 = c_0_in;

%PFR
[t, y_pfr] = ode23s('pfr', [0 tau_total], c_0);
conversion_A_pfr = (y_pfr(1,1) - y_pfr(end,1))/y_pfr(1,1);
Selectivity_B_pfr = y_pfr(end,2)/(y_pfr(1,1) - y_pfr(end,1));
yield_B_pfr = conversion_A_pfr * Selectivity_B_pfr;


%% CSTRs in series
n_cstr = 5;

% Space Time for each reactor
tau_total = tau_total/n_cstr;

% CSTR1
%Initial Concentration
c_0 = c_0_in;
y_cstr_frac1 = fsolve('cstr', c_0);

% CSTR2
%Initial Concentration
c_0 = y_cstr_frac1;
y_cstr_frac2 = fsolve('cstr', c_0);

% CSTR3
%Initial Concentration
c_0 = y_cstr_frac2;
y_cstr_frac3 = fsolve('cstr', c_0);


% CSTR4
%Initial Concentration
c_0 = y_cstr_frac3;
y_cstr_frac4 = fsolve('cstr', c_0);

% CSTR5
%Initial Concentration
c_0 = y_cstr_frac4;
y_cstr_frac5 = fsolve('cstr', c_0);

conversion_A_series = (c_0_in(1)- y_cstr_frac5(1))/c_0_in(1);
Selectivity_B_series = y_cstr_frac5(2)/(c_0_in(1)- y_cstr_frac5(1));
yield_B_series = conversion_A_series*Selectivity_B_series;

fprintf('The yield of B is : %f ', yield_B_cstr, yield_B_pfr, yield_B_series)


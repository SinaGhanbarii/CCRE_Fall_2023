% 
% - CCRE - PROF. MATTEO MAESTRI - 2023-2024 - December 11
%

clear all, clc 
% First Guess With Differential Method!!!!
global kr tau cAin cBin

Pressure=101325; %Pa
R=8.314; %J/mol/K
Q_298K_1atm=0.08; %Nm3/s
% d_r=0.02; %m reactor diameter
% L_r=0.06; %m reactor length
% A_r=pi*d_r^2/4; %m2
V_r=7.7E-6; %m3

Temperature=493; %K
Q=Q_298K_1atm*(Temperature/298); %m3/s
tau=V_r/Q; %s

% fitting procedure
firstGuess=[0.97 1]'; % kr n

cAin = [0.75
1
1.25
1.5
1.75
];

cAout = [0.728
0.97
1.214
1.456
1.699
];

cBin = [0
    0
    0
    0
    0
];

parameters=lsqcurvefit(@modelfunction,firstGuess,cAin,cAout);
figure(1)
plot(cAin,cAout,'ko',cAin,modelfunction(parameters,cAin),'b-')
legend('Data','Fitted exponential')
title('Data and Fitted Curve')

%% function
function modelfunction=modelfunction(parameters,t)
global kr n tau cAin cBin

    kr=parameters(1);
    n=parameters(2);

    cA_out = zeros(length(cAin),1);

    for i = 1:length(cAin)
        [t,y]=ode23s(@PFR,[0 tau],[cAin(i),cBin(i)]);
        cA_out(i) = y(end,1);
    end
    
    modelfunction=cA_out;

end

function PFR=PFR(t,y)
global kr n

cA = y(1);
rate = kr*cA^n;

PFR(1)=-rate;
PFR(2)=rate;

PFR=PFR';
end
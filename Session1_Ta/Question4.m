close, clear, clc
global alpha beta tau_tot Cin_tracer time y C

time = [0,0.5,1,1.5,2,2.5,3,3.5,4,5]; %min
C = [0, 0.43449,0.66865,0.81223,0.881173,0.936730,0.962299,0.978663,0.982717,0.995999]; %kmol/m3

k = 5; % m3/kmol/min
Vtot = 0.7; % l
Qtot = 0.7; % l/min
Cain = 10; % kmol/m3
tau_tot = Vtot/Qtot; % min

F = 1 - exp(-time/tau_tot);

figure
plot(time,C, 'o')
xlabel('time(min)')
ylabel('Concentration (m3/mol)')

hold on
plot(time, F)

Cin_tracer = 1;

%fitting procedure
firstGuess = [0.4 0.98]'; %alpha beta
parameters = lsqcurvefit(@modelfunction, firstGuess, time, C, [0 0], [1 1]);

figure 
plot(time,C,'ko')
hold on
plot(time, modelfunction(parameters, time),'b-')
legend('data','fitted exponential')



%% Model Function
function modelfunction = modelfunction(parameters, t)
global alpha beta tau_tot Cin_tracer time y C
    alpha = parameters(1);
    beta = parameters(2);

    C1 = ones(length(time),1);
    Cout = ones(length(time),1);

    [t,y] = ode23(@tracerEq, [0:0.5:5],0);

    for i=1:length(time)
        for j=1:length(t)
            if t(j) == time(i)
                C1(j) = y(j);
                break;
            end
        end
    end
Cout = alpha*Cin_tracer + (1-alpha)*C1;
modelfunction = Cout;
end
%%
function tracerEq = tracerEq(t,c)
    global alpha beta tau_tot Cin_tracer 
    tracerEq = 1/tau_tot * (1-alpha)/(1-beta)*(Cin_tracer-c);
end
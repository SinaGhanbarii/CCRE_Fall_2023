% --CCRE (Prof. Matteo Maestri) - 2023/2024 - December 11, 2023  
%
% A -> B 

clear all, clc

global cA time

cA = [10
9.60309042821846
9.22191873507677
8.85588187271495
8.50434802011683
8.16673008255723
7.84252759033502
7.53123069999733
7.23221663403565
6.94501892363005
6.66922957282831
6.40444058567818
6.15021733242314
5.90602013602527
5.67148382705090
5.44627251416652
5.23005030603866
5.02243938752058
4.82301168738078
4.63148221207409
4.44757443192589
4.27101181726157
4.10146216593646
3.93859776494936
3.78219021352828
3.63201185274513
3.48783438616180
3.34936527194651
3.21636183631793
3.08863665422805
2.96600230062884
2.84826682655802
2.73518147853473
2.62656500616827
2.52226310780698
2.42212148179917
2.32597575374622
2.23362182099008
2.14492157861465
2.05974799450782
1.97797403655750
1.89945681731384
1.82403444277350
1.75159913006633
1.68204629730790
1.61527136261378
1.55115625943775
1.48958306132682
1.43045727473735];

time = [0
5
10
15
20
25
30
35
40
45
50
55
60
65
70
75
80
85
90
95
100
105
110
115
120
125
130
135
140
145
150
155
160
165
170
175
180
185
190
195
200
205
210
215
220
225
230
235
240];

%integral method
firstGuess=[0.0001 2]';  

parameters=lsqcurvefit(@modelfunction,firstGuess,time,cA);
figure(1)
plot(time,cA,'ko',time,modelfunction(parameters,time),'b-')
legend('Data','Fitted exponential')
title('Data and Fitted Curve')

%% function
function modelfunction=modelfunction(parameters,t)
global k n cA time

    k=parameters(1);
    n=parameters(2);

    tspan = linspace(0,time(end),length(cA));
    [t,y]=ode23s(@batchReactor,tspan,[cA(1),0]);

    modelfunction=y(:,1);

end

function batchReactor=batchReactor(t,c)
global k n  

NS = 2; 
nu = [-1,1];
rate = k*c(1)^n;

for i = 1:NS 
   batchReactor(i) = nu(i)*rate;
end

batchReactor=batchReactor';
end































clear
close all

F=10000;
T_20_01=3.1805;
    
T_20_025= 1.0502;
T_20_05=0.4211;
T_20_075=0.2442;
T_20_1=0.1541;
T_20_2 = 0.0556;
T_50_01=3.1569;
T_50_025=1.0560;
T_50_05= 0.4126;
T_50_075=0.2331;
T_50_1=0.1552;
T_50_2= 0.0555;

   





X=size(F);
H = zeros(X);
M = [];
N = [];
n = [];

  options.N = 20;
 %options.n = 20;
 options.n=50;
 options.lambda = 1;
  options.r = 2;
 options.p = 2;

 % Generate random numbers from a half-normal distribution
rng('default'); 
%
%
%M = median(rand(20,20),2); 
% M = median(generateRandomCH(20,20,1.5),2);
% M = median(generateRandomLF(20,20,4),2);
%M = median(wblrnd(1,1.4,[50,20]),2);
%M=median(gamrnd(2,1,[50,20]),2);

for i= 1:F
  M = median(2*abs(random('Normal', 0, 2, [50, 20])),2);
[testStat,result] = TestStat_MedianExponential(M,options);
testStat  = testStat > T_50_2;
H(i) = testStat;
end 


disp(sum(H))

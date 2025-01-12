clear
close all

F=10000;
T_01= 1.9962;
T_025 = 0.5641;
T_05 = 0.2156;
T_075 = 0.1147;
T_1 = 0.0746;
T_2 = 0.0278;
X=size(F);
H = zeros(X);
S_01 = 1.3326;
S_025 = 0.3659;
S_05 = 0.1342;
S_075 = 0.0746;
S_1 = 0.0510;
S_2 = 0.0175;

options.a = 0;
options.b = 1;
options.N = 15;
options.r = 0.1;
options.p = 2;

%pd = makedist('Normal','mu',0.25,'sigma',sqrt(0.5));
%t = truncate(pd,0,1);
%
% M = median(random(t,15,15),2);

%M=median(betarnd(0.5,0.5,[10,9]),2);
for i= 1:F
M = median(generate_random_numbers_S3(2,15,15),2);
[testStat,result] = TestStat_MedianUniform(M,options);
testStat  = testStat > S_01;
H(i) = testStat;
end 


disp(sum(H))


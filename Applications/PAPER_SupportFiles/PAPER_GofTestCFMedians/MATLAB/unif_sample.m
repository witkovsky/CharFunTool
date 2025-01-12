clear 
close all
smile = [10.4, 19.6, 18.8, 13.9, 17.8, 16.8, 21.6, 17.9, 12.5, 11.1, 4.9, 12.8, 14.8, 22.8, 20.0, 15.9, 16.3, 13.4, 17.1, 14.5, 19.0, 22.8, 1.3, 0.7, 8.9, 11.9, 10.9, 7.3, 5.9, 3.7, 17.9, 19.2, 9.8, 5.8, 6.9, 2.6, 5.8, 21.7, 11.8, 3.4, 2.1, 4.5, 6.3, 10.7, 8.9, 9.4, 9.4, 7.6, 10.0, 3.3, 6.7, 7.8, 11.6, 13.8, 18.6];
% Inicijalizacija promenljivih
putnici = [1,12,4,10,4,14,11,7,11,4,13,2,4,6,3,10,0,12,6,9,10,5,13,4,10,14,12,11,6,10,11,0,11,13,2];

bstrap = [];
rezult=[];
testStat=[];
a = min(putnici);
b = max(putnici);
d= length(putnici);

options.a = min(putnici);
options.b = max(putnici);
options.N = d;
options.r = 0.1;
options.p = 2;


% Broj iteracija


n = 10000;
% Simulacija bootstrap uzoraka i testiranje


R = rand(n,options.N)*(options.b - options.a) + options.a; 
Medians = median(R,2);
Medians(n+1) = median(putnici);
testStat= TestStat_MedianUniform(Medians,options);




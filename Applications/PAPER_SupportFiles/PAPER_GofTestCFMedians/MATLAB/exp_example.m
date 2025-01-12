%set=[20, 106, 14, 78, 94, 20, 21, 136, 56, 232, 89, 33, 181, 424, 14, 430, 155, 205, 117, 253, 86, 260, 213, 58, 276, 263, 246, 341, 1105, 50, 136];
set2 = [56, 83, 104, 116, 244, 305, 429, 452, 453, 503, 552, 614, 661, 673, 683, 685, 753, 763, ...
        806, 834, 838, 862, 897, 904, 981, 1007, 1008, 1049, 1060, 1107, 1125, 1141, 1153, 1154, 1193, ...
        1201, 1253, 1313, 1329, 1347, 1454, 1464, 1490, 1491, 1532, 1549, 1568, 1574, 1586, 1599, ...
        1608, 1723, 1769, 1795, 1927, 1957, 2005, 2010, 2016, 2022, 2037, 2065, 2096, 2139, 2150, ...
        2156, 2160, 2190, 2210, 2220, 2248, 2285, 2325, 2337, 2351, 2437, 2454, 2546, 2565, 2584, ...
        2624, 2675, 2701, 2755, 2877, 2879, 2922, 2986, 3092, 3160, 3185, 3191, 3439, 3617, 3685, ...
        3756, 3826, 3995, 4007, 4159, 4300, 4487, 5074, 5579, 5623, 6869, 7739];
set_novo =[5.950, 119.077, 366.074, 155.848, 30.534, 20.615, 15.135, 3.590, 103.713, 120.859];
clear option

bstrap = [];
rezult=[];
testStat=[];
set3=[];
d = length(set_novo);
lambda=[];
N=d;
r=0.1;
p=2;

lam = 1/mean(set2);
set3 = lam * set_novo;
lambda = 1;
options.lambda = lambda;
options.N = d;
options.r = 0.1;
options.p = 2;
options.Upp = 2;

% Broj iteracija
n=10000;
R = exprnd(1/options.lambda,n,options.N);
Medians = median(R,2);
Medians(n+1) = median(set3);
testStat= TestStat_MedianExponential(Medians,options);





    



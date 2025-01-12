F=10000;
n = 30;
  N = 19;
  %M = median(rand(n,N)*2-1,2);
  M=[];
 options.a = 0;
 options.b = 1;
  options.N = N;
  alpha=0.05;
  %options.alpha   = [0.1 0.05];
  options.r= 0.25;
  %options.r = [0.1 0.15 0.25 0.4 0.5 0.6 0.75 1 2];
  options.p = 2;
 options.Upp = 5;
  [T_crit,M, options, F, alpha, n] = Median_cit(M,options,F,alpha,n);
   disp(T_crit)
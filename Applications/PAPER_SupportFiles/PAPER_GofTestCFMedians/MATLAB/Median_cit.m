function [Tcrit,M, options, F, alpha,n] =  Median_cit(M, options,F, alpha,n)
%% T_crit  
%  T_crit calculates the critical values of the GOF T statistic. 
%  
% SYNTAX:
%  [Tcrit,RND,N,a,alpha,n,tmax,funID] =  T_crit(RND,N,a,alpha,n,tmax,funID)
%
% INPUTS:
%  RND    - anonymous function which generates data from the standard bivariate
%           logistic distribution of given sample size n. If empty,
%           the function RND is loaded from the accompanying RND_2DLogistic.mat
%           file. 
%  N      - Number of simulations. Default value is N = 1000.
%  a      - vector of the weight function parameters. Default value is a =
%           2. 
%  alpha  - vector of significance levels. Default value is alpha = 0.05. 
%  n      - vector of sample sizes. Default value is n = 100. 
%  tmax   - integration limit. We set t1min = -tmax; t1max = tmax; 
%           t2min = -tmax; t2max = tmax. Default value is tmax = 5.
%  funID  - identifier of the integrand function: If funID=1 (default)
%           the integrand function is based on the ratio of the Empirical
%           CF and the null CF (Bivariate Logistic CF) suggested by Andjela
%           Mijanovic and Bozidar Popovic. If funID=2 the integrand
%           function is based on the difference of the Empirical CF and the
%           null CF (Bivariate Logistic CF). Otherwise,  the integrand
%           function is based on the difference of the Empirical CF and the
%           Bivariate Normal CF.
%
% EXAMPLE 1
% %  cf = @(t) cf2D_Logistic(t);
% %  clear options;
% %  options.isInterp = true;
% %  options.isPlot = false;
% %  result = cf2Dist2D_VW(cf,[],options);
% %  RND = result.RND;
%  RND   = [];
%  N     = 100;
%  a     = [2 3];
%  alpha = [0.1 0.05 0.01];
%  n     = [10 20];
%  [Tcrit,RND,N,a,alpha,n,tmax,funID] = T_crit(RND,N,a,alpha,n);
%  disp(Tcrit)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 03-Feb-2023 12:17:10

%% ALGORITHM
% [Tcrit,RND,N,a,alpha,n,tmax,funID] =  T_crit(RND,N,a,alpha,n,tmax,funID)

%% CHECK THE INPUT PARAMETERS
%% Algorithm

r      = options.r(:);
alpha  = alpha(:);
N      = options.N(:);
na     = length(r);
nalpha = length(alpha);
nn     = length(N);
Tcrit  =  zeros(na,nalpha,nn);
TT =  zeros(F,1);
for i = 1:na
    for j = 1:nn
        for kk = 1:F
            disp([i,j,kk])
            TT(kk)  = TestStat_MedianUniform( median(rand(n,N(j))*2-1,2),options);
        end
        T_sort = sort(TT);
        for k = 1:nalpha
            Tcrit(i,k,j) = T_sort(ceil(F*(1-alpha(k))));
        end
    end
end

end
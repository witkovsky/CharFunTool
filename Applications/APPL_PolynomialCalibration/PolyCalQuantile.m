function [q,cf] = PolyCalQuantile(x,prob,betaEstimate,LX,KY,...
    cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,order,options)
% PolyCalQuantile is an auxiliary function for the algorithm PolyCal which
% computes the quantile q of the conditional Y-distribution, for X = x, at
% the probability level specified by prob.   
%
% SYNTAX
% [q,cf] = PolyCalQuantile(x,prob,betaEstimate,LX,KY,...
%    cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,order,options)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-Sep-2020 16:46:23
%% Algorithm  

w = ones(1,order+1);
for i = 2:(order+1)
    w(i) = w(i-1).*x;
end
% CF of the predicted Y given x
cf = PolyCalCF(w*LX,w*KY,cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,w*betaEstimate);

% Quatiles of Y calculated for given prob
[~,~,~,q] = cf2DistGP(cf,[],prob,options);
end
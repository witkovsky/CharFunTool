function [xLow,xUpp] = PolyCalInterval(yLow,yUpp,betaEstimate,LX,KY,...
    xMin,xMax,probLow,probUpp,order,cfXA,cfXB,cfXB0,...
    cfYA,cfYB,cfYB0,options)
% PolyCalInterval is an auxiliary function for the algorithm PolyCal which
% calculates the calibration interval [xmin,xmax] by inverting the
% calibration function. In fact, the values xmin, xmax are solutions to the
% following optimization problems: 
%   xLow = argmin (yLow - quantile(x,probUpp))^2,
%   xUpp = argmin (yUpp - quantile(x,probLow))^2,
% where the function quantileY(x,prob) computes the quantile q of the
% conditional Y-distribution, given X = x, at the probability level
% specified by prob.  
%
% SYNTAX
%  [xLow,xUpp] = PolyCalInterval(yLow,yUpp,betaEstimate,LX,KY,...
%     xMin,xMax,probLow,probUpp,order,cfXA,cfXB,cfXB0,...
%     cfYA,cfYB,cfYB0,options)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-Sep-2020 16:46:23
%% Algorithm 
% Check ymin and ymax. 
% Valid values are from the support [yMin,yMax]

yMin = PolyCalQuantile(xMin,probLow,betaEstimate,LX,KY,...
    cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,order,options);
yMax = PolyCalQuantile(xMax,probUpp,betaEstimate,LX,KY,...
    cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,order,options);

if yLow < yMin || yUpp >yMax
    error('ymin or ymax out of range')
end

% Set the optimization function
xfun = @(x,y,pr) (y - PolyCalQuantile(x,prob,betaEstimate,LX,KY,...
    cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,order,options))^2;

% Calculate the corresponding interval limits [xmin, xmax]
xLow = fminsearch(@(x)xfun(x,yLow,probUpp),(xMax-xMin)/2);
xUpp = fminsearch(@(x)xfun(x,yUpp,probLow),(xMax-xMin)/2);

end
function [xMean,xStd] = cf2Moments(cf,tolDiff)
%cf2Moments Auxiliary function. Estimates the MEAN and the STD from the
%   characteristic function by using its numerical differentiation at 0.
%
% SYNTAX:
%  [xMean,xStd] = cf2Moments(cf,tolDiff)
%
% INPUT:
%  cf      - function handle of the characteristic function (CF), 
%  tolDiff - tolerance for numerical differentiation. By default tolDiff =
%            1e-4.
%
% EXAMPLE1 (Calculate MEAN and STD for distribution specified by CF)
%  cf = @(t) exp(-t.^2/2);
%  [xMean,xStd] = cf2Moments(cf)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-Mar-2021 11:43:31

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);

if nargin < 2, tolDiff = []; end

if isempty(tolDiff)
    tolDiff = 1e-4;
end

%%

cft   = cf(tolDiff*(1:4));
cftRe = real(cft);
cftIm = imag(cft);
xMean = (8*cftIm(1)/5 - 2*cftIm(2)/5 + 8*cftIm(3)/105 ...
        - 2*cftIm(4)/280) / tolDiff;
xM2   = (205/72 - 16*cftRe(1)/5 + 2*cftRe(2)/5 - 16*cftRe(3)/315 ...
        + 2*cftRe(4)/560) / tolDiff^2;
xStd  = sqrt(xM2 - xMean^2);

end
function [stdCF,meanCF] = PolyCalCF2Std(cf,tolDiff)
% PolyCalCF2Std is an auxiliary function for the algorithm PolyCal which
% calculates the STD and MEAN from the cell vector of characteristic
% functions. Here we assume that cf is a function handel (or cf is
% n-dimensional cell of function handles) of well defined characteristic
% function(s) of probability distribution(s) which have at leat two moments
% - the mean and the variance.
% 
% SYNTAX:
% [stdCF,meanCF] = PolyCalCF2Std(cf)
%
% EXAMPLE
% cf = {@(t)cf_Normal(t), ...
%       @(t)cf_Normal(2*t), ...
%       @(t)cf_Normal(3*t) .* cf_RectangularSymmetric(t), ...
%       @(t)cf_Student(3*t,5) .* cf_RectangularSymmetric(t), ...
%       @(t)cf_Student(5*t,5)};
% [stdCF,meanCF] = PolyCalCF2Std(cf)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-Sep-2020 16:46:23

%% Algorithm
narginchk(1, 2);
if nargin < 2, tolDiff = []; end

if isempty(tolDiff)
    tolDiff = 1e-2;
end

n = length(cf);
if n > 1
    stdCF = zeros(n,1);
    meanCF = zeros(n,1);
    m2CF = zeros(n,1);
    for i = 1:n
        cft     = cf{i}(tolDiff*(1:4));
        cftRe   = real(cft);
        cftIm   = imag(cft);
        meanCF(i)  = (8*cftIm(1)/5 - 2*cftIm(2)/5 + 8*cftIm(3)/105 ...
            - 2*cftIm(4)/280) / tolDiff;
        m2CF(i)    = (205/72 - 16*cftRe(1)/5 + 2*cftRe(2)/5 ...
            - 16*cftRe(3)/315 + 2*cftRe(4)/560) / tolDiff^2;
        stdCF(i)   = sqrt(m2CF(i) - meanCF(i)^2);
    end
elseif n==1
    cft     = cf(tolDiff*(1:4));
    cftRe   = real(cft);
    cftIm   = imag(cft);
    meanCF  = (8*cftIm(1)/5 - 2*cftIm(2)/5 + 8*cftIm(3)/105 ...
        - 2*cftIm(4)/280) / tolDiff;
    m2CF    = (205/72 - 16*cftRe(1)/5 + 2*cftRe(2)/5 ...
        - 16*cftRe(3)/315 + 2*cftRe(4)/560) / tolDiff^2;
    stdCF   = sqrt(m2CF - meanCF^2);
end
end
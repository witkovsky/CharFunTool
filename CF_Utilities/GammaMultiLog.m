function f = GammaMultiLog(z,p,funmode)
% GammaMultiLog evaluates logarithm of the multivariate gamma function of
% order p, the argument z may be complex and of any size. 
%
% Multivariate gamma function of order p, say gamma_p(z), is defined as
%  gamma_p(z) = pi^(p*(p-1)/4) * prod_{j=1}^p gamma(z-(j-1)/2).
% Hence, the logarithm of multivariate gamma is given by
%  log(gamma_p(z))=(p*(p-1)/4)*log(pi) + sum_{j=1}^p log(gamma(z-(j-1)/2))
%
% SYNTAX
%   f = GammaMultiLog(z,p,funmode)
%
% INPUTS
%  a       - complex argument, of any size (vector, matrix, array) 
%  p       - order of the multivariate gamma, if empty, default value is
%            p = 1.
%  funmode - Function mode: If funmode=0 f = log[gamma_p(z)). If funmode=1
%            f = gamma_p(z).
%
% EXAMPLE 1:
%  t = linspace(-10,10,101)';
%  z = 1i*t;
%  p = 10;
%  f = GammaMultiLog(z,p)
%
% EXAMPLE 2:
%  t = linspace(-10,10,101)';
%  z = 1i*t;
%  p = 10;
%  funmode = 1;
%  f = GammaMultiLog(z,p,funmode)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 18-Sep-2019 16:26:09

%% FUNCTION
%  f = GammaMultiLog(z,p,funmode)

%% CHECK THE INPUT PARAMETERS

if nargin < 3
    funmode = [];
end

if nargin < 2
    p = [];
end

if isempty(p)
    p = 1;
end

if isempty(funmode)
    funmode = 0;
end

sz = size(z);
z  = z(:);

% This was corrected on September 18, 2019 (from f = (p*(p-1)/2)*log(pi))
f = (p*(p-1)/4)*log(pi);

for j = 1:p
    f = f + GammaLog(z - (j-1)/2);
end

if funmode==1
    f = exp(f);
end

f = reshape(f,sz);
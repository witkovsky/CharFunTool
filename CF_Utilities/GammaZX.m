function g = GammaZX(z,funmode)
% GammaZX  
%  GAMMA function valid in the entire complex plane, the argument z
%  may be complex and of any size. 
%
% GammaLX(z) is calculated as exp(GammaLog(z)), where GammaLog calculates
% the natural Log of the Gamma function by the excellent Lanczos series
% approximation for the complex log(Gamma) function.
%
% SYNTAX
%   g = GammaZX(z,funmode)
%
% INPUTS
%  z       - complex argument, of any size (vector, matrix, array) 
%  funmode - function mode indicator. If funmod = 1 (the default value)
%            GammaZX(z,funmode=1) = exp(GammaLog(z)). If funmod = 0,
%            GammaZX(z,funmode=0) = GammaLog(z).

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 25-Sep-2019 13:56:14

%% FUNCTION CALL
%  g = GammaZX(z,funmode)

%% ALGORITHM

if nargin < 2
    funmode = [];
end

if isempty(funmode)
    funmode = 1;
end

if funmode == 0
    g = GammaLog(z);
else
    g = exp(GammaLog(z));
end
end
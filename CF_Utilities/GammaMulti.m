function f = GammaMulti(z,p)
% GammaMulti evaluates the multivariate gamma function of
% order p, the argument z may be complex and of any size. 
%
% Multivariate gamma function of order p, say gamma_p(z), is defined as
%  gamma_p(z) = pi^(p*(p-1)/2) * prod_{j=1}^p gamma(z-(j-1)/2).
%
% SYNTAX
%   f = GammaMulti(z,p)
%
% INPUTS
%  a       - complex argument, of any size (vector, matrix, array) 
%  p       - order of the multivariate gamma, if empty, default value is
%            p = 1.
%
% EXAMPLE 1:
%  t = linspace(-10,10,101)';
%  z = 1i*t;
%  p = 10;
%  f = GammaMulti(z,p)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Aug-2018 19:45:37

%% FUNCTION
%  f = GammaMulti(z,p)

%% ALGORITHM

f = GammaMultiLog(z,p,1);

end
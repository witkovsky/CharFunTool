function cf = cfS_Student(t,df,mu,sigma,coef,niid)
%cfS_Student
%  Characteristic function of the STUDENT's t-distribution with df > 0
%  degrees of freedom.
%
%  cfS_Student is an ALIAS of the more general function cf_Student, used to
%  evaluate the characteristic function of a linear combination of
%  independent (location-scale) STUDENT's t-distributed random variables.
%
%  The characteristic function of the STUDENT's t-distribution with df
%  degrees of freedom is defined by
%   cf(t) = cfS_Student(t,df) 
%         = besselk(df/2,abs(t)*sqrt(df),1) * exp(-abs(t)*sqrt(df)) * ...
%          (sqrt(df)*abs(t))^(df/2) / 2^(df/2-1)/gamma(df/2);
%
% SYNTAX:
%  cf = cfS_Student(t,df)
%  cf = cfS_Student(t,df,mu,sigma,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df    - the degrees of freedom, df > 0. If empty, the default value is
%          df = 1. 
%  mu    - vector of location parameters, mu in Real. If empty, default
%          value is mu = 0. 
%  sigma - vector of scale parameters, sigma_i > 0. If empty, default value
%          is sigma = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient n, such that Z = Y + ... + Y is
%          sum of niid random variables Y, where each Y = sum_{i=1}^N
%          coef_i * X_i is independently and identically distributed random
%          variable. If empty, default value is niid = 1. 
% 
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Student%27s_t-distribution.
%
% EXAMPLE 1:
% % CF of the Student t-distribution with df = 3
%   df = 3;
%   t = linspace(-5,5,501);
%   cf = cfS_Student(t,df);
%   figure; plot(t,cf),grid
%   title('CF of the Student t-distribution with df = 3)')
%
% EXAMPLE 2:
% % PDF/CDF of the Student t-distribution with df = 3
%   df = 3;
%   cf = @(t) cfS_Student(t,df);
%   x = linspace(-8,8,101);
%   prob = [0.9 0.95 0.975 0.99];
%   clear options
%   options.N = 2^12;
%   options.SixSigmaRule = 30;
%   result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
%   WITKOVSKY, V.: On the exact computation of the density and of the
%   quantiles of linear combinations of t and F random variables. Journal
%   of Statistical Planning and Inference 94 (2001), 1–13.

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Jun-2017 12:08:24
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
narginchk(1, 6);
if nargin < 6,  niid = []; end
if nargin < 5,  coef = []; end
if nargin < 4, sigma = []; end
if nargin < 3,    mu = []; end
if nargin < 2,    df = []; end

cf = cf_Student(t,df,mu,sigma,coef,niid);

end
function cf = cfS_Student(t,df)
%cfS_Student
%  Characteristic function of the STUDENT's t-distribution with df > 0
%  degrees of freedom.
%
%  cfS_Student is a version resp. ALIAS of the more the more general
%  function cf_Student, used to evaluate the characteristic function of a
%  linear combination of independent  STUDENT's t-distributed random
%  variables.
%
%  The characteristic function of the STUDENT's t-distribution with df
%  degrees of freedom is defined by
%   cf(t) = cfS_Student(t,df) 
%         = besselk(df/2,abs(t)*sqrt(df),1) * exp(-abs(t)*sqrt(df)) * ...
%          (sqrt(df)*abs(t))^(df/2) / 2^(df/2-1)/gamma(df/2);
%
% SYNTAX:
%  cf = cfS_Student(t,df);
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  df     - the degrees of freedom, df > 0. If empty, the default value is
%           df = 1. 
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
%   clear options
%   options.SixSigmaRule = 8;
%   result = cf2DistGP(cf,x,[],options)
%
% REFERENCES:
%   WITKOVSKY, V.: On the exact computation of the density and of the
%   quantiles of linear combinations of t and F random variables. Journal
%   of Statistical Planning and Inference 94 (2001), 1–13.

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Jun-2017 12:08:24

%% ALGORITHM
narginchk(1, 2);
if nargin < 2,  df = []; end
if isempty(df), df = 1; end

cf = cf_Student(t,df);
end
function f = Hypergeom1F1MatApprox(a,b,X)
%Hypergeom1F1MatApprox Computes the approximation of the confluent
%  hypergeometric function 1F1(a,b,X) of a matrix argument, defined
%  for the complex parameters a and b, with Re(a) > (p-1)/2 and  Re(b-a) >
%  (p-1)/2, and a REAL symmetric (p x p)-matrix argument X.
%
%  In fact, 1F1(a,b,X) depends only on the eigenvalues of X, so X could be
%  specified as a (p x p)-diagonal matrix or a p-dimensional vector of
%  eigenvalues of the original matrix X, say x.
%
%  Based on heuristic arguments (not formally proved yet), the value of the
%  confluent hypergeometric function 1F1(a,b,X) of a matrix argument is
%  calculated as  
%   1F1(a;b;X) ~ 1F1(a;b;x(1)) * ... * 1F1(a;b;x(p)),
%  where 1F1(a;b;x(1)) is the scalar value confluent hypergeometric
%  function 1F1(a,b,x(i)) with [x(1),...,x(p)] = eig(X). 
%
%  Here the confluent hypergeometric function 1F1(a;b;z) is evaluated for
%  the vector parameters a and b and the scalar argument z by using the
%  simple (4-step) series expansion. 
%
% SYNTAX:
%  f = Hypergeom1F1MatApprox(a,b,X)
%
% INPUTS:
%  a      - complex vector of parameters of the hypergeometric function
%           1F1(a;b;X),
%  b      - complex vector of parameters  of the hypergeometric function
%           1F1(a;b;X),
%  X      - real symmetric (p x p)-matrix argument (alternatively can be
%           specified as a (p x p)-diagonal matrix or a p-vector of the
%           eigenvalues of X), i.e. x = eig(X).
%
% OUTPUT:
%  f     - (Approximate) value of the confluent hypergeometric function
%           1F1(a;b;X), of a matrix argument X. 
%
% EXAMPLE 1:
%  a = 3;
%  b = 5;
%  X = [1,2];
%  f = Hypergeom1F1MatApprox(a,b,X)
%
% EXAMPLE 2:
% % PDF/CDF of minus log Wilks Lambda RV (p=10, n=30, q=5) from its CF
% % Here, cf_LogRV_WilksLambdaNC id based on using Hypergeom1F1MatApprox
%   p      = 10;
%   n      = 30;
%   q      = 5;
%   Delta  = [1 2 3 10 50]; % nonzero eigenvalues of non-centrality matrix
%   coef   = -1;
%   cf     = @(t) cf_LogRV_WilksLambdaNC(t,p,n,q,Delta,coef);
%   prob   = [0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result  = cf2DistGP(cf,[],prob,options)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 19-Jul-2018 16:11:57

%% FUNCTION CALL
%  f = Hypergeom1F1MatApprox(a,b,X)

%% CHECK THE INPUT PARAMETERS

sza = size(a);
szb = size(b);
sz = [max([sza(1),szb(1)]),max([sza(2),szb(2)])];

[errorcode,a,b] = distchck(2,a(:),b(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% ALGORITHM

[p1,p2] = size(X);
if p1 == p2
    x = eig(X);
elseif min(p1,p2) == 1
    x = X;
else
    error(message('InputSizeMismatch'));
end

p = length(x);
f = 1;
for i = 1:p
    f = f .* Hypergeom1F1SeriesExp(a,b,x(i));
end

f = reshape(f, sz);
end
%% Function Hypergeom1F1SeriesExp
function [f,isConv,loops,n,tol] = Hypergeom1F1SeriesExp(a, b, x, n, tol)
%Hypergeom1F1SeriesExp Computes the confluent hypergeometric function
%  1F1(a;b;z) also known as the Kummer's function M(a,b,z), for the
%  vector parameters a and b and the scalar argument z by using the simple
%  (4-step) series expansion.
%
%  For more details on confluent hypergeometric function 1F1(a;b;z) or the
%  Kummer's (confluent hypergeometric) function M(a, b, z) see WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Confluent_hypergeometric_function.
%
% SYNTAX
%   f = Hypergeom1F1SeriesExp(a, b, x, n)
%   [f,isConv,loops,n,tol] = Hypergeom1F1SeriesExp(a, b, x, n, tol)
%
% INPUTS
%  a      - vector parameter a,
%  b      - vector parameter b,
%  x      - scalar argument x,
%  n      - maximum number of terms used in the series expansion. If empty,
%           default value of n = 500,
%  tol    - tolerance for stopping rule. If empty, default value f tol =
%           1e-14.
%
% OUTPUTS
%  f      - calculated 1F1(a,b,x) of the same dimension as a snd b,
%  isConv - flag indicator for convergence of the series, isConv = true if
%           loops < n
%  loops  - total number of terms used in the series expansion,
%  n      - used maximum number of terms in the series expansion,
%  tol    - used tolerance for stopping rule.
%
% EXAMPLE 1
%  t = linspace(-20,20,101);
%  a = 1i*t;
%  b = 2 - 3i*t;
%  x = 3.14;
%  f = Hypergeom1F1SeriesExp(a, b, x)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 19-Jul-2018 15:15:51

%% FUNCTION
%  [f,isConv,loops,n,tol] = Hypergeom1F1SeriesExp(a, b, x, n, tol)

%% CHECK THE INPUT PARAMETERS

narginchk(3, 5);

if nargin < 4, n   = []; end
if nargin < 5, tol = []; end

if isempty(n)
    n = 500;
end

if isempty(tol)
    tol = 1e-14;
end

sza = size(a);
szb = size(b);

if prod(sza) >= prod(szb)
    sz = sza;
else
    sz = szb;
end

[errorcode,a,b] = distchck(2,a(:),b(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end

f  = 1;
r1 = 1;
for j  = 1 : 4 : n
    r1 = r1 .* (a + j - 1) ./ ( j * ( b + j - 1)) * x;
    r2 = r1 .* (a + j) ./ ( (j + 1) * ( b + j)) * x;
    r3 = r2 .* (a + j + 1) ./ ( (j + 2) * ( b + j + 1)) * x;
    r4 = r3 .* (a + j + 2) ./ ( (j + 3) * ( b + j + 2)) * x;
    rg = r1 + r2 + r3 + r4;
    f  = f + rg;
    if( max(abs( rg ./ f )) < tol )
        break;
    end
    r1 = r4;
end

loops = j;
isConv = false;
if loops < n
    isConv = true;
end

f = reshape(f,sz);
end
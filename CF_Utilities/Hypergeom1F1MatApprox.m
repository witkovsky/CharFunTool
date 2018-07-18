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
%  Here, based on heuristic arguments, the approximate value of 1F1(a;b;X)
%  is calculated as 
%   1F1(a;b;X) ~ 1F1(a;b;x(1)) * ... * 1F1(a;b;x(p)), 
%  where x = eig(X).
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
%  f     - Approximae confluent hypergeometric function 1F1(a;b;X).  
%
% EXAMPLE 1:
%  a = 3;
%  b = 5;
%  X = [1,2];
%  f = Hypergeom1F1MatApprox(a,b,X)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 18-Jul-2018 22:35:51

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
    f = f .* Hypergeom1F1(a,b,x(i));
end

f    = reshape(f, sz);
end
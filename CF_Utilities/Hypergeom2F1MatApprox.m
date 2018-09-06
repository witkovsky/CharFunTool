function f = Hypergeom2F1MatApprox(a,b,c,X)
%Hypergeom2F1MatApprox Computes approximation of the Gauss hypergeometric
%  function 2F1(a,b,c,X) of a (p x p)-matrix argument X. Hypergeom2F1Mat is
%  defined for the complex  parameters a, b, and c with Re(a) > (p-1)/2 and
%  Re(c-a) > (p-1)/2, and  a real symmetric matrix argument X, with Re(X) <
%  I. 
%
%  For more details and definition of the hypergeometric functions
%  with matrix argument see, e.g., Koev and Edelman (2006) or Muirhead
%  (2009). 
% 
%  In fact, 2F1(a,b,c,X) depends only on the eigenvalues of X, so X could
%  be specified as a (p x p)-diagonal matrix or a p-dimensional vector of
%  eigenvalues of the original matrix X, say x.
%
%  Based on heuristic arguments (not formally proved), the approximate
%  value of the  Gauss hypergeometric function 2F1(a,b;c;X) of a matrix
%  argument is calculated here as  
%   2F1(a,b,c,X) ~ 2F1(a,b,c,x(1)) * ... * 2F1(a,b,c,x(p)),
%  where 2F1(a,b,c,x(i)) is the scalar value of the Gauss hypergeometric
%  function 2F1(a,b,c,x(i)) with [x(1),...,x(p)] = eig(X). 
%
% SYNTAX:
%  f = Hypergeom2F1MatApprox(a,b,c,X)
%
% INPUTS:
%  a      - complex vector of parameters of the hypergeometric function
%           2F1(a,b,c,X),
%  b      - complex vector of parameters  of the hypergeometric function
%           2F1(a,b,c,X)),
%  c      - complex vector of parameters  of the hypergeometric function
%           2F1(a,b,c,X)),
%  X      - real symmetric (p x p)-matrix argument (alternatively can be
%           specified as a (p x p)-diagonal matrix or a p-vector of the
%           eigenvalues of X), i.e. x = eig(X).
%
% OUTPUT:
%  f     - (Approximate) value of the Gauss hypergeometric function
%           2F1(a,b,c,X)) of a matrix argument X. 
%
% EXAMPLE 1:
%  a = 3;
%  b = 5;
%  c = 4;
%  X = [0.5,0.75];
%  f = Hypergeom2F1MatApprox(a,b,c,X)
%
% EXAMPLE 2
%  t = 1i*linspace(-5,5,11)';
%  a = 3*t;
%  b = 2.5*t;
%  c = 1.5*t;
%  X = [0.5,0.75];
%  f = Hypergeom2F1MatApprox(a,b,c,X)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 05-Sep-2018 22:55:20

%% FUNCTION CALL
%  f = Hypergeom2F1MatApprox(a,b,c,X)

%% CHECK THE INPUT PARAMETERS


sza = size(a);
szb = size(b);
szc = size(c);

[errorcode,a,b,c] = distchck(3,a(:),b(:),c(:));
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
%    f = f .* Hypergeom2F1(a,b,c,x(i));
    f = f .* HypergeompFqSeries([a,b],c,x(i));
end

if max(sza)>1
    f = reshape(f,sza);
    return
elseif max(szb)>1
    f = reshape(f,szb);
    return
elseif max(szc)>1
    f = reshape(f,szc);
end
end
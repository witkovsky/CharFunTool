function f = Hypergeom2F1Mat(a,b,c,X,MAX)
% Hypergeom2F1Mat - The Gauss hypergeometric function 2F1(a,b;c;X) of a
%  (p x p)-matrix argument X. Hypergeom2F1Mat is defined for the complex
%  parameters a, b, and c with Re(a) > (p-1)/2 and  Re(c-a) > (p-1)/2, and
%  a real symmetric matrix argument X, with Re(X) < I. 
%
%  For more details and definition of the hypergeometric functions
%  with matrix argument see, e.g., Koev and Edelman (2006) or Muirhead
%  (2009). 
%
% SYNTAX:
%  f = Hypergeom2F1Mat(a,b,c,X,MAX)
%
% INPUTS:
%  a      - complex vector of parameters of the hypergeometric function
%           2F1(a,b;c;X),
%  b      - complex vector of parameters  of the hypergeometric function
%           2F1(a,b;c;X),
%  c      - complex vector of parameters  of the hypergeometric function
%           2F1(a,b;c;X), 
%  X      - real symmetric (p x p)-matrix argument (alternatively can be
%           specified as a (p x p)-diagonal matrix or a p-vector of the
%           eigenvalues of X).
%  MAX    - maximum number of partitions, |kappa|<=MAX, default value is
%           MAX = 20,
%
% OUTPUTS:
%  f     - hypergeometric sum, 2F1(a,b;c;X)
%
% EXAMPLE:
%  a   = 3;
%  b   = 2.5;
%  c   = 1.5;
%  X   = [1,2,3]/5;
%  MAX = 50;
%  f   = Hypergeom2F1Mat(a,b,c,X,MAX)
%
% REFERENCES:
% [1] Koev, P. and Edelman, A., 2006. The efficient evaluation of the
%     hypergeometric function of a matrix argument. Mathematics of
%     Computation, 75(254), 833-846.
% [2] Muirhead RJ. Aspects of multivariate statistical theory. John Wiley &
%     Sons; 2009 Sep 25.
% [3] Butler RW, Wood AT. Laplace approximations for hypergeometric
%     functions with matrix argument. The Annals of Statistics.
%     2002;30(4):1155-77.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 22-Aug-2018 12:23:08

%% FUNCTION CALL
% f = Hypergeom2F1Mat(a,b,c,X,MAX)

%% CHECK THE INPUT PARAMETERS
if nargin < 5
    MAX = 20; 
end

sza = size(a);
szb = size(b);
szc = size(c);

[errorcode,a,b,c] = distchck(3,a(:),b(:),c(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end



%% ALGORITHM
f = HypergeompFqMat([a,b],c,X,[],2,MAX);

if max(sza)>1
    f = reshape(f,sza);
    return
elseif max(szb)>1
    f = reshape(f,szb);
    return
elseif max(szc)>1
    f = reshape(f,szc);
end
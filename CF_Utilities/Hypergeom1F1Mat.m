function f = Hypergeom1F1Mat(a,b,X,MAX)
% Hypergeom1F1Mat - The confluent hypergeometric function 1F1(a;b;X) of a
%  (p x p)-matrix argument X. Hypergeom1F1Mat is defined for the complex
%  parameters a and b, with Re(a) > (p-1)/2 and  Re(b-a) > (p-1)/2, and a
%  symmetric matrix argument X, specified as a vector of its eigenvalues.
%
%  For more details and definition of the hypergeometric functions
%  with matrix argument see, e.g., Koev and Edelman (2006) or Muirhead
%  (2009). 
%
% SYNTAX:
%  [s,ss] = Hypergeom1F1Mat(a,b,X,MAX)
%
% INPUTS:
%  a      - complex vector of parameters of the hypergeometric function
%           1F1^alpha(a;b;X),  
%  b      - complex vector of parameters  of the hypergeometric function
%           1F1^alpha(a;b;X), 
%  X      - matrix argument (specified as vector of its eigenvalues),
%  MAX    - maximum number of partitions, |kappa|<=MAX, default value is
%           MAX = 20,
%
% OUTPUTS:
%  f     - hypergeometric sum, 1F1(a;b;X)
%
% EXAMPLE 1
%  a   = 3;
%  b   = 5;
%  X   = [1,2,3];
%  MAX = 10;
%  f   = Hypergeom1F1Mat(a,b,X,MAX)
%
% EXAMPLE 2
%  a   = [3 4 5]';
%  b   = [5+1i 5+2i 5+3i]';
%  X   = [1,2,3];
%  MAX = 50;
%  f   = Hypergeom1F1Mat(a,b,X,MAX)
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
% Ver.: 25-Sep-2019 15:55:45

%% ALGORITHM

if nargin < 4
    MAX = 20; 
end

sza = size(a);
szb = size(b);

[errorcode,a,b] = distchck(2,a(:),b(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end

try
    % MEX-FUNCTION CALL
    f = HypergeompFqMat_mex(a,b,X,[],2,MAX,[]);
catch
    % M-FUNCTION CALL
    f = HypergeompFqMat(a,b,X,[],2,MAX,[]);
end

if max(sza)>1
    f = reshape(f,sza);
    return
elseif max(szb)>1
    f = reshape(f,szb);
end
end
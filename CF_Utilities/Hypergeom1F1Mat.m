function f = Hypergeom1F1Mat(a,b,X,MAX)
% Hypergeom1F1Mat - The confluent hypergeometric function 1F1(a;b;X) of a
%  (p x p)-matrix argument X. Hypergeom1F1Mat is defined for the complex
%  parameters a and b, with Re(a) > (p-1)/2 and  Re(b-a) > (p-1)/2, and a
%  REAL symmetric matrix argument X.  
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
%  X      - real symmetric (p x p)-matrix argument (alternatively can be
%           specified as a (p x p)-diagonal matrix or a p-vector of the
%           eigenvalues of X).
%  MAX    - maximum number of partitions, |kappa|<=MAX, default value is
%           MAX = 20,
%
% OUTPUTS:
%  f     - hypergeometric sum, 1F1(a;b;X)
%
% EXAMPLE:
%  a   = 3;
%  b   = 5;
%  X   = [1,2,3];
%  MAX = 10;
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
% Ver.: 25-Oct-2017 14:56:37

%% FUNCTION CALL

if nargin < 4
    MAX = 20; 
end

f = HypergeompFqMat(a,b,X,[],2,MAX,[]);

end
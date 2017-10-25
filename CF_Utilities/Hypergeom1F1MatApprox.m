function f = Hypergeom1F1MatApprox(a,b,X)
%Hypergeom1F1MatApprox Computes the Laplace approximation of the confluent
%  hypergeometric function 1F1(a,b,X) of a (p x p)-matrix argument, defined
%  for the complex parameters a and b, with Re(a) > (p-1)/2 and  Re(b-a) >
%  (p-1)/2, and a REAL symmetric matrix argument X. 
%
%  In fact, 1F1(a,b,X) depends only on the eigenvalues of X, so X could be
%  specified as a (p x p)-diagonal matrix or a p-dimensional vector of
%  eigenvalues of the original matrix X. 
%
% SYNTAX:
%  f = Hypergeom1F1MatApprox(a,b,X)
%
% INPUTS:
%  a      - complex vector of parameters of the hypergeometric function
%           1F1^alpha(a;b;X),  
%  b      - complex vector of parameters  of the hypergeometric function
%           1F1^alpha(a;b;X), 
%  X      - real symmetric (p x p)-matrix argument (alternatively can be
%           specified as a (p x p)-diagonal matrix or a p-vector of the
%           eigenvalues of X).
%
% OUTPUT:
%  f     - Laplace approximation of the confluent hypergeometric function
%          1F1(a;b;X), calibrated such that 1F1(a;b;0) = 1.  
%
% EXAMPLE 1:
%  a = 3;
%  b = 5;
%  X = [1,2];
%  f = Hypergeom1F1MatApprox(a,b,X)
%
% REFERENCES:
% [1] Butler RW, Wood AT. Laplace approximations for hypergeometric
%     functions with matrix argument. The Annals of Statistics.
%     2002;30(4):1155-77.
% [2] Butler RW, Wood AT. Approximation of power in multivariate analysis.
%     Statistics and Computing. 2005 Oct 1;15(4):281-7. 

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Oct-2017 22:15:55

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
f = p*(b-(p+1)/4) .* log(b);
y = zeros(length(a),p);
for i = 1:p
    D = (x(i)-b).^2 + 4*a.*x(i);
    y = 2*a ./ (b - x(i) + sqrt(D));
    f = f +  a.*log(y./a) + (b-a).*log((1-y)./(b-a)) + x(i)*y;
    y(:,i) = y;
end

r = 1;
for i = 1:p
    for j = i:p
        r = r .* ( y(:,i).*y(:,j)./a + (1-y(:,i)).*(1-y(:,j))./(b-a) );
    end
end

f = reshape( exp(f-log(r)/2), sz);

end
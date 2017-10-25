function f = Hypergeom2F1MatApprox(a,b,c,X)
%Hypergeom2F1MatApprox Computes the Laplace approximation of the Gauss
%  hypergeometric functionfunction 2F1(a,b;c;X) of a (p x p)-matrix
%  argument X. Hypergeom2F1MatApprox is defined for the complex parameters
%  a, b, and c with Re(a) > (p-1)/2 and  Re(c-a) > (p-1)/2, and a symmetric
%  matrix argument X, with Re(X) < I. 
%
%  In fact, 2F1(a,b;c;X) depends only on the eigenvalues of X, so X could
%  be specified as a (p x p)-diagonal matrix or a p-dimensional vector of 
%  eigenvalues of the original matrix X. 
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
%
% OUTPUT:
%  f     - Laplace approximation of the Gauss hypergeometric function
%          2F1(a,b;c;X), calibrated such that 2F1(a,b;c;0) = 1.    
%
% EXAMPLE:
%  a   = 10;
%  b   = 2.5;
%  c   = 5;
%  X   = [1,2,3]/5;
%  f   = Hypergeom2F1MatApprox(a,b,c,X)
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
%  f = Hypergeom2F1MatApprox(a,b,c,X)

%% CHECK THE INPUT PARAMETERS
sza = size(a);
szb = size(b);
szc = size(c);
sz = [max([sza(1),szb(1),szc(1)]),max([sza(2),szb(2),szc(2)])];

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
f = p*(c-(p+1)/4) .* log(c);
y = zeros(length(a),p);
S = zeros(length(a),p);
for i = 1:p
    tau = x(i)*(b-a) - c;   
    D = tau.^2 - 4*a.*x(i).*(c-b);
    y = 2*a ./ (sqrt(D) - tau);
    f = f +  a.*log(y./a) + (c-a).*log((1-y)./(c-a)) - b.*log(1-x(i)*y);
    S(:,i) = x(i).*y.*(1-y) ./ (1-x(i).*y);
    y(:,i) = y;
end

r = 1;
for i = 1:p
    for j = i:p
        r = r .* ( y(:,i).*y(:,j)./a + (1-y(:,i)).*(1-y(:,j))./(c-a) ...
            - b.*S(:,i).*S(:,j)./(a.*(c-a)) );
    end
end

f = reshape( exp(f-log(r)/2), sz);

end
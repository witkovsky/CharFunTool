function [roots,warning,err] = FindRoots(fun,A,B,n,isplot,tol,maxiter)
%FINDROOTS estimates the real roots (zeros) of a real (oscillatory)
% function FUN on the interval [A,B], by using adaptive nth-order (n=2^k)
% Chebyshev polynomial approximation of the function FUN.
%
% This code was adapted from the code published in Day & Romero (2005):
% Roots Of Polynomials Expressed In Terms Of Orthogonal Polynomials. SIAM
% Journal on Numerical Analysis,  43, 1969 - 1987.
%
% SYNTAX:
% roots = FindRoots(fun,A,B)
% roots = FindRoots(fun,A,B,n,isplot)
% [roots,warning,err] = FindRoots(fun,A,B,n,isplot,tol)
%
% INPUTS:
%  fun    - function handle, e.g. fun = @(x)sin(x)
%  A,B    - lower and upper limit of the interval [A,B]. Default values are
%           A = -1, and B = 1;
%  n      - polynomial order for approximation of the function fun, on each 
%           (automatically selected) subinterval of [A,B]. n should be 
%           power of 2. Default value is n = 2^5.  
%           Increase the value of n if FindRoots is unable to find all 
%           roots of fun over the interval [A,B].
%  isplot - logical flag. If isplot = true, FindRoots plots the graph of
%           the function together with depicted locations of its roots.
%  tol    - tolerance for stoping the rootFinder. If empty, default value
%           is tol = 1e-16.
%
% EXAMPLE 1:
% fun = @(t) sin(t.^3 + t.^2 + t)
% A   = 0; 
% B   = 5;
% roots = FindRoots(fun,A,B)
%
% EXAMPLE 2:
% fun = @(t)  exp(-.3*t) .* sin(10*t) .* cos(2*t.^2)
% A   = 5; 
% B   = 15;
% roots = FindRoots(fun,A,B)
%
% EXAMPLE 3:
% x   = 3;
% nu  = 1;
% cf_chi2 = @(t) (1 - 2i * t) .^(-nu/2);
% fun = @(t) min(4,imag(exp(-1i*t*x).*cf_chi2(t))./t)
% A   = 0.2; 
% B   = 10;
% n   = 2^7;
% roots = FindRoots(fun,A,B,n)
%
% EXAMPLE 4:
% nu  = 3;
% fun = @(t)  sin(0.5+5*t) .* (besselj(nu,t) - besselk(nu,t))
% A   = 150; 
% B   = 200;
% roots = FindRoots(fun,A,B)
%
% REFERENCES:
%   Day, D., Romero, L. (2005): Roots Of Polynomials Expressed In Terms Of
%   Orthogonal Polynomials. SIAM Journal on Numerical Analysis, 43,
%   1969 - 1987.
%
% ACKNOWLEDGEMENS:
% Chebfun inspired this file. For more details see http://www.chebfun.org/

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jul-2017 10:06:48

%% FUNCTION
%  [roots,warning,err] = FindRoots(fun,A,B,n,isplot)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 7);
if nargin < 7, maxiter = []; end
if nargin < 6, tol = []; end
if nargin < 5, isplot = []; end
if nargin < 4, n = []; end
if nargin < 3, B = []; end
if nargin < 2, A = []; end

%% SET THE DEFAULT VALUES (input parameters)

if isempty(maxiter)
    maxiter = 100;
end

if isempty(tol)
    tol = 1e-16;
end

if isempty(isplot)
    isplot = true;
end

if isempty(n)
    n = 2^5;
end

if isempty(B)
    B = 1;
end

if isempty(A)
    A = -1;
end

interval = [A;B];

if abs(B-A) > 3*pi
    divisionrule = [-0.5;0;0.5];
    interval =  GetSubs(interval,divisionrule);
else
    divisionrule = 0;
end

%% ALGORITHM
[cgl, M] = setupChebyshev(n);

roots = [];
warning = 0;
iter = 0;
while iter < maxiter
    iter = iter +1;
    [r,interval,err,isWarning] = rootFinder(cgl,M,fun,interval,n,tol);
    roots = cat(1,roots,r);
    if isempty(interval)
        break
    else
        interval =  GetSubs(interval,divisionrule);
    end
    warning = warning + isWarning;
end

if iter == maxiter
    warning = warning + 1;
end

roots = unique(roots);

%% Plot the function together with the roots
if isplot
    N = 1000;
    t = linspace(A,B,N);
    plot(t,fun(t)),grid
    hold on
    plot(roots,fun(roots),'or')
    hold off
end
end % End of function FindRoots
%% Function setupChebyshev
function [cgl,M] = setupChebyshev(n)
% setupChebyshev evaluates the Chebyshev-Gauss-Legendre points and the
% Chebyshev Transformation Matrix
%
% Adapted from the MATLAB code published in Day & Romero (2005): Roots of
% polynomials expressed in terms of orthogonal polynomials. SIAM Journal on
% Numerical Analysis,  43, 1969 - 1987.

% Chebyshev-Gauss-Legendre points
ind = (0:n);
cgl = cos(ind * pi/n);

% Chebyshev Transformation Matrix
M        = cos(ind'* ind * pi/n);
M(1,:)   = M(1,:)/2;
M(:,1)   = M(:,1)/2;
M(n+1,:) = M(n+1,:)/2;
M(:,n+1) = M(:,n+1)/2;
M        = M * (2/n);
end % End of function setupChebyshev
%% Function evalCheb
function V = evalCheb(n,z)
%EVALCHEB evaluates the required Vandermonde matrix
%
% Input: vector of points, z, and the polynomial degree n.
% Output: Vandermonde matrix, m by degree_max + 1, V(j+1,k)=T_j(z_k)
%
% Adapted from the MATLAB code published in Day & Romero (2005): Roots of
% polynomials expressed in terms of orthogonal polynomials. SIAM Journal on
% Numerical Analysis,  43, 1969 - 1987.

m = size(z,1);
if m*n >= 0
    V(:,1) = ones(m,1);
    if n >= 1
        V(:,2) = z;
        if n >= 2
            index = find( log(abs(z)) >= 100/n );
            si = size(index,1);
            if si > 0
                z(index) = ones(si,1)*exp(100/n);
            end
            for i=2:n
                V(:,i+1) = V(:,i).*(2*z) - V(:,i-1);
            end
        end
    end
else
    V = [];
end
end % End of function evalCheb
%%
function [roots,intervals,err,isWarning] = rootFinder(cgl,M,fun,intervals,n,tol)
%ROOTFINDER estimates the roots of fun over intervals by using nth order
% Chebyshev polynomial approximation of fun
%
% Adapted from the MATLAB code published in Day & Romero (2005): Roots of
% polynomials expressed in terms of orthogonal polynomials. SIAM Journal on
% Numerical Analysis,  43, 1969 - 1987.

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 12-Aug-2013 16:18:42

roots = [];

nint = size(intervals,2);
err  = false(nint,1);
isWarning = 0;
for i = 1:nint
    range = (intervals(2,i) - intervals(1,i))/2;
    shift = (intervals(1,i) + intervals(2,i))/2;
    x = cgl*range + shift;
    f = fun(x);
    ExpansionCoeff = f * M;    
    if abs(ExpansionCoeff(n+1)) < tol
        isWarning = 1;
        % warning('The leading expansion coefficient vanishes');        
    else
        ExpansionCoeff = ExpansionCoeff/(-2*ExpansionCoeff(n+1));
        H      = diag(ones(n-1, 1)/2, 1) + diag(ones(n-1, 1)/2, -1);
        H(1,2) = 1;
        C      = H;
        C(n,:) = C(n,:) + ExpansionCoeff(1:n);    
        Eigenvalues = eig(C);
        Vandermonde = evalCheb(n,Eigenvalues);
        Vcolsums    = sum(abs(Vandermonde'));
        tube_index  = find((abs(imag(Eigenvalues))<.2) ...
                      & (abs(real(Eigenvalues))< 2));
        Solutions   = Eigenvalues(tube_index);
        Vcolsums    = Vcolsums(tube_index);
        cond_max    = min(2^(n/2), 10^6);
        condEigs_index =  Vcolsums < cond_max ;
        Solutions   = sort(Solutions(condEigs_index));
        r = range * Solutions + shift;
        if isreal(r)
            r = r(r>=intervals(1,i) & r<=intervals(2,i));
            roots = cat(1,roots,r);
        else
            err(i) = true;
        end
    end
end
intervals = intervals(:,err);
end % End of function rootFinder
%% FUNCTION GETSUBS (Sub-division of the integration intervals)
function SUBS = GetSubs(SUBS,XK)
%GETSUBS Sub-division of the intervals for adaptive root finding
%
% Example
% interval = [0;1];
% divisionrule = [-.25; -.5; 0; .5; .75] % selected division points from [-1;1]
% intervals = GetSubs(interval,divisionrule)
% intervals = GetSubs(intervals,divisionrule)

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 20-May-2013 01:28:00

NK = size(XK,1);
SUBIND = 1:NK;
M = 0.5*(SUBS(2,:)-SUBS(1,:));
C = 0.5*(SUBS(2,:)+SUBS(1,:));
Z = XK*M + ones(NK,1)*C;

L = [SUBS(1,:); Z(SUBIND,:)];
U = [Z(SUBIND,:); SUBS(2,:)];
SUBS = [reshape(L, 1, []); reshape(U, 1, [])];
end % End of function GetSubs
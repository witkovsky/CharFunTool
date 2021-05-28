function [x, w] = ChebPoints(n,domain)
% ChebPoints evaluates n Chebyshev points (of the 2nd kind) in the given
%  interval (domain = [a,b]), together with the row vector of the weights
%  for Clenshaw-Curtis quadrature (over the given interval [a,b]).
% 
%  The ChebPoints algorithm is based on the ChebFun chebpts algorithm. For
%  more details see http://www.chebfun.org/.  
%  
% SYNTAX:
%   [x, w] = ChebPoints(n,domain)
%
% EXAMPLE1 (Barycentric interpolant of the Sine function on (-pi,pi))
%   domain = [-pi,pi];
%   x      = ChebPoints(32,domain);
%   f      = sin(x);
%   xNew   = linspace(domain(1),domain(2),201)';
%   fBary  = InterpBarycentric(x,f,xNew);
%   fCheb  = InterpChebValues(f,xNew,domain)
%   fTrue  = sin(xNew);
%   plot(x,f,xNew,fBary,xNew,fCheb,'.')
%   disp([fTrue fBary fCheb])
%
% EXAMPLE2 (Integral of the Sine function on the interval (0,pi))
%   domain = [0,pi];
%   [x,w] = ChebPoints(32,domain);
%   f     = sin(x);
%   Itrue = 2;
%   Icalc = w * f;
%   disp([Itrue Icalc])

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-May-2021 14:28:24
% Revisions: 24-Jul-2017 10:06:48
%
% Based on the algorithm CHEBPTS of ChebFun.
% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%% FUNCTION
%  [x, w] = ChebPoints(n,domain)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);
if nargin < 2, domain = []; end

if isempty(domain)
    domain = [-1,1];
end

%% ALGORITHM

if n == 0
    x = [];
elseif n == 1
    x = 0;
else
    m = n - 1;
    x = sin(pi*(-m:2:m)/(2*m)).';
end

if domain(1) ~= -1 || domain(2) ~= 1
    x = domain(2)*(x + 1)/2 + domain(1)*(1 - x)/2;
end

if nargout > 1
    if n == 0
        w = [];
    elseif n == 1
        w = 2;
    else
        v = 2./[1, 1-(2:2:(n-1)).^2];  % Integrals of ChebPoly T_k (even)
        v = [v, v(floor(n/2):-1:2)];   % Mirror for DCT via FFT
        w = ifft(v);                   % Interior weights
        w([1,n]) = w(1)/2;             % Boundary weights
    end   
    if domain(1) ~= -1 || domain(2) ~= 1
        w = (diff(domain)/2)*w;
    end
end

end
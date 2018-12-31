function [fun,cdfFun,pdfFun] = IntegrandFun_BTAV(phi,x,cf,funtype,M)
%IntegrandFun_BTAV 
%  Auxiliary function calculates the integrand functions for computing the
%  CDF/PDF of the NON-NEGATIVE DISTRIBUTION specified by its characteristic
%  function, by using the Bromwich-Talbot-Abate-Valko (BTAV) inversion
%  method (originally suggested as numerical inversion for the Laplace
%  transform function). In particular, the Bromwich contour is properly
%  deformed and the integrand functions are given by
%    pdfFun = const * exp(s0.*x) .* cf(1i*s0) .* s1;
%    cdfFun = pdf ./ s0;
%  where const = 2*M/5 (default value is M = 10), cf is an anonymous
%  characteristic function of nonnegative distribution which is well
%  defined for complex valued arguments, and 
%    s0 = const * phi .* (cot(phi)+1i) ./ x
%    s1 = const * 1i * (1 + 1i*(phi + (phi.*cot(phi)-1).*cot(phi))) ./ x
%  for -pi <= phi <= pi.
%  Then the CDF/PDF can be numerically evaluated at specified values x by
%  integrating the integrand functions over the interval (-pi,pi) - by
%  using any quadrature rule (e.g. the simple trapezoidal or the more
%  advanced adaptive Gauss-Kronrod quadrature rule):
%    CDF = real(integral(cdfFun,-pi,pi,'ArrayValued',true))/(2*pi)
%    PDF = real(integral(pdfFun,-pi,pi,'ArrayValued',true))/(2*pi)
%  or alternatively
%    CDF = integral(@(phi)real(cdfFun(phi)),0,pi,'ArrayValued',true))/(pi)
%    PDF = integral(@(phi)real(pdfFun(phi)),0,pi,'ArrayValued',true))/(pi)
%  For more details see Talbo (1979) and Abate & Valko (2004).
%
% SYNTAX:
%  [fun,cdfFun,pdfFun] = IntegrandFun_BTAV(phi,x,cf,funtype,M)
%
% EXAMPLE 1
% % CDF/PDF of chi-square distribution with 1 degree of freedom, df = 1
% % By using integral - Matlab adaptive Gauss-Kronrod quadrature
%  x      = [2.705543454095416 6.634896601021214 28.373987362798132];
%  cf     = @(t)(1-2i*t).^(-1/2);
%  cdfFun = @(phi) IntegrandFun_BTAV(phi,x,cf,'cdf');
%  pdfFun = @(phi) IntegrandFun_BTAV(phi,x,cf,'pdf');
%  CDF    = real(integral(cdfFun,-pi,pi,'ArrayValued',true))/(2*pi);
%  PDF    = real(integral(pdfFun,-pi,pi,'ArrayValued',true))/(2*pi);
%
% EXAMPLE 2 
% % CDF/PDF of chi-square distribution with 1 degree of freedom, df = 1
% % By using simple trapezoidal quadrature
%  x      = [2.705543454095416 6.634896601021214 28.373987362798132];
%  cf     = @(t)(1-2i*t).^(-1/2);
%  N      = 100;
%  phi    = linspace(-pi,pi,N+1)';
%  [~,cdfFun,pdfFun] = IntegrandFun_BTAV(phi(2:end-1),x,cf);
%  CDF    = real(sum(cdfFun))/N;
%  PDF    = real(sum(pdfFun))/N;
%
% EXAMPLE 3 
% % CDF/PDF of chi-square distribution with 1 degree of freedom, df = 1
% % By using efficient trapezoidal quadrature on half interval
%  x      = [2.705543454095416 6.634896601021214 28.373987362798132];
%  cf     = @(t)(1-2i*t).^(-1/2);
%  N      = 50;
%  phi    = linspace(0,pi,N+1)';
%  [~,cdfFun,pdfFun] = IntegrandFun_BTAV(phi,x,cf);
%  CDF    = (cdfFun(1,:)/2 + sum(real(cdfFun(2:end-1,:))))/N;
%  PDF    = (pdfFun(1,:)/2 + sum(real(pdfFun(2:end-1,:))))/N;
%
% REFERENCES:
% [1] Talbot, A., 1979. The accurate numerical inversion of Laplace
%     transforms. IMA Journal of Applied Mathematics, 23(1), pp.97-120.
% [2] Abate, J. and Valkó, P.P., 2004. Multi-precision Laplace transform
%     inversion. International Journal for Numerical Methods in
%     Engineering, 60(5), pp.979-993.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 25-Dec-2018 12:42:35

%% ALGORITHM
%[fun,cdfFun,pdfFun] = IntegrandFun_BTAV(phi,x,cf,funtype,M)

%% CHECK THE INPUT PARAMETERS
narginchk(3, 5);

if nargin < 5, M = []; end
if nargin < 4, funtype = []; end

if isempty(M)
    M = 10;
end

if isempty(funtype)
    funtype = 'cdf';
end

const = 2*M/5;
phi   = phi(:);
x     = x(:)';
ct    = cot(phi);
s0    = phi .* (ct+1i);
s0(phi==0) = 1;
s0    = const * s0 ./ x;
s1    = 1 + 1i*(phi+(phi.*ct-1).*ct);
s1(phi==0) = 1;
s1    = s1 ./ x;

pdfFun = const * exp(s0.*x) .* cf(1i*s0) .* s1;
cdfFun = pdfFun ./ s0;

pdfFun(abs(phi)==pi) = 0;
cdfFun(abs(phi)==pi) = 0;

switch lower(funtype)
    case 'cdf'
        fun = cdfFun;
    case 'pdf'
        fun = pdfFun;
end
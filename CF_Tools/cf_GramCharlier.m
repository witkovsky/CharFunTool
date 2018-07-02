function [cf,an,dn,kn] = cf_GramCharlier(t,kn,kn0,cf0,options)
%% cf_GramCharlier evaluates the characteristic function of the approximate
%  probability distribution specified by its moments (cumulants) based on
%  the generalized Gram-Charlier expansion (Type A series expansion).
%
%  The usual form of the Gram-Charlier expansion is an expansion about a
%  normal distribution with common mean mu and standard error sigma (i.e.,
%  common cumulants k1 and k2). A generalization of the Gram-Charlier
%  expansion allows to express CF of one distribution in terms of another.
%  The coefficients of this general expansion are explicitly obtained from
%  the given cumulants. For more details see Berberan-Santos (2007) and
%  notes in the REMARK.
%
%  In particular, let kn denote the cumulants of the specified distribution
%  of interest (for n = 1,2,...) and let cf0(t) is the characteristic
%  function of the reference distribution (e.g. standard normal
%  distribution specified by its cumulants kn0 for n = 1,2,...). The the
%  characteristic function of the distribution specified by the cumulants
%  kn is
%   cf(t) = sum_{n=0}^Inf (an*(1i*t)^n/n!) * cf0(t),
%  where the coeficients an (for n = 1,2,...) are uniquely derived from the
%  cumulant differences dn = kn - kn0, n = 1,2,... . For specification of
%  the coefficient, see the REMARKS section below.
%
% SYNTAX
%  cf = cf_GramCharlier(t,kn)
%  cf = cf_GramCharlier(t,kn,kn0,cf0,options)
%
% INPUTS:
%  t       - vector or array of real values, where the CF is evaluated.
%  kn      - vector of cumulants or moments of the specified distribution
%            of interest (this should be specified by the options). By
%            default kn is assumed to be vector of cumulants, kn for n =
%            1,...,nMax.
%  kn0     - vector of cumulants or moments of the reference distribution
%            (this should be specified by the options). By default kn0 is
%            assumed to be vector of cumulants, kn0 for n = 1,...,nMax. If
%            empty, kn0 = [ 0,1,0,0,...]', i.e. a vector of cumulants of
%            the standard normal N(0,1) distribution.
%            CHECK that the stated cumulants kn0 are  cumulants of the
%            specified reference characteristic function cf0(t) (MUST BE!).
%  cf0     - function handle of the reference characteristic function. If
%            empty, default value the reference characteristic function is
%            the standard normal distribution, cf0(t)= @(t)cfS_Gaussian(t).
%  options - structure with the following default parameters:
%            options.isMoment = false.
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Edgeworth_series.
%
% EXAMPLE1 (CF of the approaximate distribution specified by kn = [0 1 0 1])
%  kn  = [0 1 0 1];
%  cf  = @(t) cf_GramCharlier(t,kn);
%  t   = linspace(-10,10,501);
%  plot(t, real(cf(t)),t,imag(cf(t)));
%  xlabel('t')
%  ylabel('CF')
%  title('CF of the distribution specified by the cummulants')
%
% EXAMPLE2 (PCF/CDF/QF of the approaximate distribution with kn = [0 1 0 1])
%  kn     = [0 1 0 1];
%  cf     = @(t) cf_GramCharlier(t,kn);
%  x      = linspace(-5,5,201);
%  prob   = [0.9 0.95 0.99];
%  options.N = 2^12;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE2 (PCF/CDF/QF of the approaximate distribution specified by its cumulants)
%  kn     = [0 1 0 1 1 1];
%  cf     = @(t) cf_GramCharlier(t,kn);
%  x      = linspace(-5,5,201);
%  prob   = [0.9 0.95 0.99];
%  options.N = 2^12;
%  result = cf2DistGP(cf,x,prob,options)
%
% REMARKS:
%  It is well known that with the exception of the degenerate Dirac
%  distribution (known also as the delta distribution) and normal
%  distribution, all distributions have an infinite number of non-zero
%  cumulants. The cumulant generating function is defined as C(t) =
%  log(cf(t)), whose series expansion gives
%   C(t) = sum_{n=1}^Inf kn *(1i*t)^n/n!,
%  where the kn are the cumulants.
%
%  Under relatively general conditions, the moments (cumulants) of a
%  distribution define the respective distribution (PDF). It is therefore
%  of interest to know how to build CF of the associated approximate
%  probability distribution specified by its moments (cumulants). An
%  obvious practical application is to obtain an approximate form from a
%  finite set of moments (cumulants).
%
%  Berberan-Santos (2007) suggested explicit formula relating a PDF with
%  its cumulants and suggested a generalization of the Gram-Charlier
%  expansion that allows to express a PDF in terms of another PDF.
%
%  The raw moments are explicitly related to the cumulants by
%   m1 = k1,
%   m2 = k1^2 + k2,
%   m3 = k1^3 + 3*k1*k2 + k3,
%   m4 = k1^4 + 6*k1^2 + 3*k2^2 + 4*k1*k3 + k4,
%   m5 = k1^5 + 10*k1^3*k2 + 15*k1*k2^2 + 10*k1^2*k3 + ...
%             + 10*k2*k3 + 5*k1*k4 + k5,
%   m6 = k1^6 + 15*k1^4*k2 + 45*k1^2*k2^2 + 15*k2^3 +  ...
%             + 20*k1^3*k3 + 60*k1*k2*k3 + 10*k3^2 + ...
%             + 15*k1^2*k4 + 15*k2*k4 + 6*k1*k5 + k6.
%  In general, the moments are defined by cumulant by the following formula
%   m{n+1} = sum_{p=0}^{n} choosenk(n,p)*m{n-p}*k{p+1}.
%
%  On the other hand, the first four cumulants are
%   k1 = m1 = mu,
%   k2 = m2 - m1^2 = sigma^2,
%   k3 = 2*m1^3 - 3*m1*m2 + m3 = gamma1*sigma^3,
%   k4 = -6*m1^4 + 12*m1^2*m2 - 3*m2^2 - 4*m1*m3 + m4 = gamma2*sigma^4,
%  where mu is the mean, sigma is the standard error, gamma1 is the
%  skewness, and gamma2 is the kurtosis. In general,
%   k{n+1} = m{n+1} - sum_{p=0}^{n-1} choosenk(n,p)*m{n-p}*k{p+1}.
%
%  Let dn = kn - kn0, for n = 1,2,..., denote the cumulant differences
%  between the specified distribution of interest and the reference
%  distribution. Then we define the associated coefficients an by
%   a0 = 1,
%   a1 = d1,
%   a2 = d1^2 + d2,
%   a3 = d1^3 + 3*d1*d2 + d3,
%   a4 = d1^4 + 6*d1^2*d2 + 3*d2^2+ 4*d1*d3 + d4,
%   a5 = d1^5 + 10*d1^3*d2 + 15*d1*d2^2 + 10*d1^2*d3 + ...
%             + 10*d2*d3 + 5*d1*d4 + d5,
%   a6 = d1^6 + 15*d1^4*d2 + 45*d1^2*d2^2 + 15*d2^3 + 20*d1^3*d3 + ...
%             + 60*d1*d2*d3 + 10*d3^2 + 15*d1^2*d4 + 15*d2*d4 + ...
%             + 6*d1*d5 + d6.
%  In general, the coeficients are defined by the cumulant differences by
%  the following formula
%   a{n+1} = sum_{p=0}^{n} choosenk(n,p)*a{n-p}*d{p+1}.
%  In particular, the coefficients of the Gram-Charlier expansion, where
%  the reference distribution is the standard normal N(0,1), and with d1 =
%  d2 = 0, are given by
%   a0 = 1,
%   a1 = 0,
%   a2 = 0,
%   a3 = d3 = k3,
%   a4 = d4 = k4,
%   a5 = d5 = k5,
%   a6 = 10*d3^2 + d6 = 10*k3^2 + k6.
%
% REFERENCES:
%  Berberan-Santos, M.N. (2007). Expressing a probability density function
%  in terms of another PDF: A generalized Gram-Charlier expansion. Journal
%  of Mathematical Chemistry, 42(3), pp.585-594.

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Jul-2018 13:57:05
%% ALGORITHM
%cf = cf_GramCharlier(t,kn,kn0,cf0,options)

%% CHECK THE INPUT PARAMETERS
narginchk(2, 5);
if nargin < 5, options = []; end
if nargin < 4, cf0 = []; end
if nargin < 3, kn0 = []; end

kn   = kn(:);
nMax = length(kn);

if ~isfield(options, 'isMoment')
    options.isMoment = false;
end

if isempty(kn0)
    kn0 = zeros(nMax,1);
    kn0(2) = 1;
    if isempty(cf0)
        cf0 = @(t) cfS_Gaussian(t);
    else
        warning('CHECK if the cumulants and the reference CF are corresponding!')
    end
end

if ~isempty(kn0) && isempty(cf0)
    kn0   = kn0(:);
    mu    = kn0(1);
    sigma = sqrt(kn0(2));
    kn0(3:end) = 0;
    cf0 = @(t) cfS_Gaussian(sigma*t).*exp(1i*t*mu);
end

if options.isMoment
    mn     = kn;
    mn0    = kn0;
    kn     = zeros(nMax,1);
    kn0    = kn;
    kn(1)  = mn(1);
    kn0(1) = mn0(1);
    for n = 2:nMax
        aux  = 0;
        aux0 = 0;
        for p = 0:n-1
            c = nchoosek(n,p) ;
            aux  = aux + c * mn(n-p) * kn(p+1);
            aux0 = aux + c * mn0(n-p) * kn0(p+1);
        end
        kn(n)   = mn(n)  - aux;
        kn0(n)  = mn0(n) - aux0;
    end
end


% Set the cummulant differences and the coefficients an
an    = zeros(nMax+1,1);
dn    = zeros(nMax+1,1);
dn(1) = 1;
dn(2:end) = kn - kn0;
an(1) = dn(1);
an(2) = dn(2);
for n = 2:nMax
    auxa = 0;
    for p = 0:n
        c = nchoosek(n,p) ;
        auxa  = auxa + c * dn(1+n-p) * an(p+1);
    end
    an(n+1)  = auxa;
end

%% Characteristic function
szt   = size(t);
t     = t(:);
nFac  = 1;
cf    = cf0(t);
kf    = log(cf);
logit = log(1i*t);
for n = 1:nMax
    nFac = nFac*n;
    cf  = cf + exp(kf + log(an(n+1)) + n*logit - log(nFac));
end
cf = reshape(cf,szt);

end
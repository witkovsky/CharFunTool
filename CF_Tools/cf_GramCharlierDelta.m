function cf = cf_GramCharlierDelta(t,delta,cfRef)
%% cf_GramCharlierDelta evaluates the characteristic function (CF) of the 
%  approximate probability distribution specified by its cumulant
%  differences, based on the generalized Gram-Charlier expansion (Type A 
%  series expansion).
%
%  The usual form of the Gram-Charlier expansion is an expansion about a
%  normal distribution with common mean mu and standard deviation sigma
%  (i.e., common cumulants kappa_1 and kappa_2). A generalization of the
%  Gram-Charlier expansion allows to express CF of one distribution in
%  terms of another. The coefficients of this general expansion are
%  explicitly obtained from the given cumulants and or moments. For more
%  details see Berberan-Santos (2007) and notes in the REMARK.
%
%  Characteristic function of the (approximate) distribution specified by
%  the cumulant differences delta_n = kappa_n - kappaRef_n for n =
%  1,...,nMax is 
%   cf(t) = exp(sum_{n=0}^nMax (delta_n*(1i*t)^n/n!)) * cfRef(t),
%  or alternatively
%   cf(t) = (sum_{n=0}^Inf (alpha_n*(1i*t)^n/n!)) * cfRef(t),
%  where delta_n = (kappa_n - kappaRef_n) and the coeficients alpha_n, for
%  n = 1,...,nMax, are uniquely derived from the cumulant differences
%  delta_n. For specification of the coefficient, see the REMARKS section
%  below. 
%
% SYNTAX
%  cf = cf_GramCharlierDelta(t,delta)
%  cf = cf_GramCharlierDelta(t,delta,cfRef)
%
% INPUTS:
%  t        - vector or array of real values, where the CF is evaluated.
%  delta    - vector of cumulant differences delta_n between the specified
%             distribution of interest and the reference distribution. 
%  cfRef    - function handle of the reference characteristic function. If
%             empty, default value the reference characteristic function is
%             the standard normal distribution, i.e. cfRef(t) = @(t)
%             cfS_Gaussian(t).
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Edgeworth_series
%
% EXAMPLE1 (CF of the approaximate distribution given with delta=[0 0 0 1])
%  delta  = [0 0 1 -1];
%  cf  = @(t) cf_GramCharlierDelta(t,delta);
%  t   = linspace(-10,10,501);
%  plot(t, real(cf(t)),t,imag(cf(t)));
%  xlabel('t')
%  ylabel('CF')
%  title('CF of the distribution specified by the cummulants')
%
% EXAMPLE2 (PDF/CDF/QF of the distribution specified with delta)
%  delta  = [0 0 1 -1];
%  cf  = @(t) cf_GramCharlierDelta(t,delta);
%  x      = linspace(-5,5,201);
%  prob   = [0.9 0.95 0.99];
%  clear options
%  options.N = 2^12;
%  options.SixSigmaRule = 10;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF/QF of the approaximate distribution specified by delta)
%  delta  = [0 0 5 1 1 1 -1 -1];
%  sigma  = 2;
%  cfRef  = @(t) cfS_Gaussian(t*sigma);
%  cf     = @(t) cf_GramCharlierDelta(t,delta,cfRef);
%  x      = linspace(-10,10,201);
%  prob   = [0.9 0.95 0.99];
%  clear options
%  options.N = 2^12;
%  options.SixSigmaRule = 10;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE4 ((CF of the approaximate distribution with ChiSquare reference)
%  delta  = [0 0 5 1 1 1 -1 -1];
%  df     = 20;
%  cfRef  = @(t) cfX_ChiSquare(t,df);
%  cf     = @(t) cf_GramCharlierDelta(t,delta,cfRef);
%  t      = linspace(-1,1,501);
%  plot(t, real(cf(t)),t,imag(cf(t)));
%  xlabel('t')
%  ylabel('CF')
%  title('CF of the distribution specified by the cummulants and cfRef')
%
% EXAMPLE5 (PDF/CDF/QF of the approaximate distribution with ChiSquare ref)
%  delta  = [0 0 5 1 1 1 -1 -1];
%  df     = 20;
%  cfRef  = @(t) cfX_ChiSquare(t,df);
%  cf     = @(t) cf_GramCharlierDelta(t,delta,cfRef);
%  x      = linspace(0,50,201);
%  prob   = [0.9 0.95 0.99];
%  clear options
%  options.N = 2^12;
%  options.SixSigmaRule = 10;
%  result = cf2DistGP(cf,x,prob,options)
%
% REMARKS:
%  It is well known that with the exception of the degenerate Dirac
%  distribution and the normal distributions, all other distributions have
%  an infinite number of non-zero cumulants. The cumulant generating
%  function is defined as C(t) = log(cf(t)), whose series expansion gives
%   C(t) = sum_{n=1}^Inf kappa_n *(1i*t)^n/n!,
%  where the kappa_n are the cumulants.
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
%  expansion that allows to express a PDF in terms of another PDF. Moreover
%  the author presented explicit formula relating the CF. 
%
%  In general, such CF (specified by a subset of cumulants) may NOT be a
%  proper characteristic function of a probability function (frequently, CF
%  is divergent for large values t). The used series is not necessarily
%  convergent, even if the standard normal is used as the reference
%  distribution. Convergence is assured for square integrable PDFs, when
%  the PDF expansion based on the derivatives of reference distribution
%  PDFs coincides with a generalized Fourier series whose basis functions
%  are orthogonal eigenfunctions of a regular Sturm-Liouville problem. For
%  more details and other references see Berberan-Santos (2007).
%
%  The raw moments are explicitly related to the cumulants (for simplicity
%  here denoted by k1, k2, k3, k4, k5, k6) by 
%   m1 = k1,
%   m2 = k1^2 + k2,
%   m3 = k1^3 + 3*k1*k2 + k3,
%   m4 = k1^4 + 6*k1^2 + 3*k2^2 + 4*k1*k3 + k4,
%   m5 = k1^5 + 10*k1^3*k2 + 15*k1*k2^2 + 10*k1^2*k3 + ...
%             + 10*k2*k3 + 5*k1*k4 + k5,
%   m6 = k1^6 + 15*k1^4*k2 + 45*k1^2*k2^2 + 15*k2^3 +  ...
%             + 20*k1^3*k3 + 60*k1*k2*k3 + 10*k3^2 + ...
%             + 15*k1^2*k4 + 15*k2*k4 + 6*k1*k5 + k6.
%  In general, the moments are defined recursively by using the known
%  cumulants by the following formula 
%   m{n+1} = sum_{p=0}^{n} choosenk(n,p)*m{n-p}*k{p+1}.
%
%  On the other hand, the first four cumulants are
%   k1 = m1 = mu,
%   k2 = m2 - m1^2 = sigma^2,
%   k3 = 2*m1^3 - 3*m1*m2 + m3 = gamma1*sigma^3,
%   k4 = -6*m1^4 + 12*m1^2*m2 - 3*m2^2 - 4*m1*m3 + m4 = gamma2*sigma^4,
%  where mu is the mean, sigma is the standard deviation, gamma1 is the
%  skewness, and gamma2 is the kurtosis. In general,
%   k{n+1} = m{n+1} - sum_{p=0}^{n-1} choosenk(n,p)*m{n-p}*k{p+1}.
%
%  Let delta_n = kappa_n - kappaRef_n denote the cumulant differences for
%  n = 1,...,nMax between the specified distribution of interest and the
%  reference distribution. The coefficients alpha_n are uniquely derived
%  from delta_n as follows (for simplicity here we use an instead of
%  alpha_n and dn instead of delta_n):
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
% [1] Berberan-Santos, M.N. (2007). Expressing a probability density
%     function in terms of another PDF: A generalized Gram-Charlier
%     expansion. Journal of Mathematical Chemistry, 42(3), pp.585-594.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 18-Sep-2018 11:56:47

%% ALGORITHM
%cf = cf_GramCharlierDelta(t,delta,cfRef)

%% CHECK THE INPUT PARAMETERS
narginchk(2, 3);
if nargin < 3, cfRef = []; end

if isempty(cfRef)
        cfRef = @(t) cfS_Gaussian(t);
end
   
% Set the coefficients alpha_n for n = 0,1,2,...
nMax     = length(delta);
alpha    = zeros(nMax+1,1);
delta        = [1;delta(:)];
alpha(1) = delta(1);
alpha(2) = delta(2);
for n = 1:(nMax-1)
    aux = 0;
    for p = 0:n
        c = nchoosek(n,p) ;
        aux  = aux + c * alpha((n-p)+1) * delta((p+1)+1);
    end
    alpha((n+1)+1)  = aux;
end

%% Characteristic function
szt   = size(t);
t     = t(:);
nFac  = 1;
cf    = cfRef(t);
kf    = log(cf);
logit = log(1i*t);
for n = 1:nMax
    nFac = nFac*n;
    cf  = cf + exp(kf + log(alpha(n+1)) + n*logit - log(nFac));
end
cf = reshape(cf,szt);

end
function cf = cfX_LogNormal(t,mu,sigma,tol)
%cfX_LogNormal(t,mu,sigma) 
% Characteristic function cf(t) of the Lognormal distribution with
% parameters mu (real) and sigma > 0, computed for real (vector) argument
% t, i.e. 
%   cf(t) = cfX_LogNormal(t,mu,sigma);
%
% NOTE: 
% This is an experimental algorithm which is not fully documented and that
% may provide unreliable results for a particular combination of parameters
% and / or may cause an unexpected failure.
%
% WIKIPEDIA:
% https://en.wikipedia.org/wiki/Log-normal_distribution
%  
% SYNTAX:
%  cf = cfX_LogNormal(t,mu,sigma)
%  cf = cfX_LogNormal(t,mu,sigma,tol)
%
% EXAMPLE1 (CF of the Lognormal distribution with mu = 0,sigma = 1)
%  mu = 0;
%  sigma = 1;
%  t = linspace(-20,20,2^10+1)';
%  cf = cfX_LogNormal(t,mu,sigma);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Lognormal distribution')
%
% EXAMPLE2 (CDF/PDF of the Lognormal distribution with mu = 0,sigma = 1)
%  mu = 0;
%  sigma = 1;
%  x = linspace(0,15,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.N = 2^13;
%  options.SixSigmaRule = 8;
%  cf = @(t) cfX_LogNormal(t,mu,sigma);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Poisson-Lognormal distribution)
%  mu = 0;
%  sigma = 1;
%  lambda = 10;
%  cfX = @(t)cfX_LogNormal(t,mu,sigma);
%  cf = @(t) cfN_Poisson(t,lambda,cfX);
%  x = linspace(0,70,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% FURTHER DETAILS:
% In probability theory, a log-normal (or lognormal) distribution is a
% continuous probability distribution of a random variable whose logarithm
% is normally distributed. The lognormal distribution is defined for x in
% (0,+inf) by its PDF/CDF/CF, as follows
%  pdf(x) = \frac{1}{x\sigma\sqrt{2\pi}}\ 
%            e^{-\frac{\left(\ln x-\mu\right)^2}{2\sigma^2}},
%  cdf(x) = \frac12 + \frac12\,\mathrm{erf}\Big[\frac{\ln x-\mu}
%            {\sqrt{2}\sigma}\Big],
%  cf(t)  = \sum_{n=0}^{\infty}\frac{(it)^n}{n!}e^{n\mu+n^2\sigma^2/2}.
% As noted, this representation is asymptotically divergent but sufficient
% for numerical purposes. 
%
% cfX_LogNormal is based on the standard integral representation of the
% characteristic function of the lognormal distribution, i.e.
%  cf(t) = Integral_0^inf exp(i*t*x) * PDF(x) dx.
% By using the half-space Fourier integral transformation we get
%  cf(t) = Integral_0^inf (i/t) * exp(-x) * PDF(i*x/t) dx.
% If we define the integrand as funCF(t,x) = (i/t) * exp(-x) * PDF(i*x/t),
% then by using a stabilizing transformation from [0,inf] to [0,1], we can
% evaluate the CF by the following (well behaved) integral:
%  cf(t) = Integral_0^1 2x/(1-x)^3 * funCF(t,(x/(1-x))^2) dx.
%
% cfX_LogNormal evaluates this integral by using the MATLAB built in
% function 'integral', with precission specified by tolerance tol (default
% value is tol = 1e-6). 
%
% REFERENCES:
% [1] WITKOVSKY, V.: On the exact computation of the density and of
%     the quantiles of linear combinations of t and F random
%     variables. Journal of Statistical Planning and Inference 94
%     (2001), 1–13.
% [2] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-Sep-2020 13:25:21
%
% Revision history:
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfX_LogNormal(t,mu,sigma,tol);

%% CHECK THE INPUT PARAMETERS

if nargin < 1
    error(message('VW:cfX_LogNormal:TooFewInputs'));
end
if nargin < 2, mu = []; end
if nargin < 3, sigma = []; end
if nargin < 4, tol = []; end

if isempty(mu), mu = 0; end
if isempty(sigma), sigma = 1; end
if isempty(tol), tol = 1e-6; end

sigma(sigma <= 0) = NaN;
reltol = tol;

%% EVALUATE THE CHARACTERISTIC FUNCTION: cfX_LogNormal(t,mu,sigma)

sz = size(t);
t  = t(:);
cf = ones(size(t));
id = t~=0;
cf(id) = integral(@(x) bsxfun(@times,funCF(mu,sigma,t(id),(x/(1-x))^2), ...
    2*x/(1-x)^3),0,1,'ArrayValued',true,'RelTol',reltol);
cf = reshape(cf,sz);

end

%% Function funCF
function f = funCF(mu,sigma,t,x)
%funCF Integrand function of the integral representation of the
%  characteristic function CF of the Lognormal distribution with (scalar) 
%  parameters mu (real) and sigma > 0 and the real (vector) argument t.
%
% SYNTAX:
%   f = funCF(mu,sigma,t,x)
%
% EXAMPLE: (Plot the integrand functions for computing the CF(t))
%  mu = 0;
%  sigma = 1;
%  t   = 1:5;
%  x   = linspace(0,1);
%  f   = funCF(mu,sigma,t,x);
%  plot(x,real(f),x,imag(f))
%
% REFERENCES:
%  WITKOVSKY V. (2016). On computing the characteristic functions of 
%  lognormal distribution and its applications. Working Paper.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 11-Apr-2016 14:06:18

%% ALGORITHM FUNCF
t  = exp(mu)*t(:);
x  = x(:)';
ot = ones(size(t));
ox = ones(size(x));
funPDF  = @(x,s) exp(-0.5*(log(x)/s).^2) ./ (sqrt(2*pi)*s*x);
if (sigma>=1/3)
    t = 1./t;
    f = (1i*t*ox) .* exp(-ot*x) .* funPDF(1i*t*x,sigma);
else
    % Set optimum limits for small and large abs(t)
    small = 7*sqrt(1/sigma);
    large = 25*sqrt(1/sigma);
    f = ot * ox;
    id = (abs(t)<=small);
    if any(id)
        f(id,:) = exp(1i*t(id)*x) .* funPDF(ot(id)*x,sigma);
    end
    id = (t>small & t<=large);
    if any(id)
        f(id,:) = exp(1i.*ot(id)*x) .* ... 
            exp(-0.5*(log((1./t(id))*x)./sigma).^2) ./ ...
            (sqrt(2*pi)*sigma*ot(id)*x);
    end
    id = (t<-small & t>=-large);
    if any(id)
        f(id,:) = exp(-1i*ot(id)*x) .* ...
            exp(-0.5*(log(-(1./t(id))*x)/sigma).^2) ./ ...
            (sqrt(2*pi)*sigma*ot(id)*x);
    end
    id = (abs(t)>large);
    if any(id)
        % TODO: Find efficient algorithm / approximation
        % From practical point of view we can set f(t,x) = 0 for all x
        f(id,:) = 0;
        % Alternatively,
        % f(id,:) = (1i*(1./t(id))*ox) .* exp(-ot(id)*x) .* ...
        %           funPDF(1i*(1./t(id))*x,sigma);
    end
end

end
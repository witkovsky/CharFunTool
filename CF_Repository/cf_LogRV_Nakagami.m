function cf = cf_LogRV_Nakagami(t,m,Omega,coef,niid)
%cf_LogRV_Nakagami 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED NAKAGAMI random variables.     
%  
%  The Nakagami distribution (or the Nakagami-m distribution) is a
%  continuous probability distribution for positive-valued random
%  variables, related to the chi-distribution, respectively. 
%
%  In particular, if X ~ Nakagami(m,Omega), with the real shape parameter m
%  >= 1/2 and the real spread parameter Omega > 0. Then, the following
%  relationship holds true: X = sqrt(Omega/2*m) * Y, where Y ~ Chi(2*m).
%
%  cf_LogRV_Nakagami evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * log(X_i), where X_i ~ Nakagami(m_i,Omega_i) are
%  inedependent Nakagami RVs with the shape parameters m_i >= 1/2 and the
%  spread parameters Omega_i > 0, for i  = 1,...,N. 
%
%  The characteristic function of log(X) with X ~ Nakagami(m,Omega) is
%   cf_LogRV_Nakagami(t) = exp(1i*t*log(sigma)) * cf_LogRV_Chi(t,df)
%  where sigma = sqrt(Omega/2*m) and cf_LogRV_Chi(t,df) denotes the
%  characteristic function of the log-transformed Chi distribution with
%  df=2*m degrees of freedom.
%
%  Hence, the characteristic function of Y is 
%   cf(t) = Prod ( cf_LogRV_Nakagami(coef_i*t,m_i,Omega_i) )
%
% SYNTAX:
%  cf = cf_LogRV_Nakagami(t,m,Omega,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  m     - vector of the shape parameters of the Nakagami random
%          variables, m >= 1/2.  If m is scalar, it is assumed that all
%          shape parameters are equal. If empty, default value is scale =
%          1/2.  
%  Omega - vector of the spread parameters of the Nakagami random
%          variables, Omega > 0.  If Omega is scalar, it is assumed that
%          all spread parameters are equal. If empty, default value is
%          Omega = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          Nakagami random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Nakagami_distribution
%  https://en.wikipedia.org/wiki/Chi_distribution
%
% NOTES (from Wikipedia)
%  The Nakagami distribution was first proposed by Nakagami in1960. It has
%  been used to model attenuation of wireless signals traversing multiple
%  paths and to study the impact of fading channels on wireless
%  communications. The PDF of the Nakagami distributioon is
%   pdf(x) = 2m^m / (Gamma(m)*Omega^m) * x^{2m-1} * exp(-{(m/Omega)*x^2),
%  for all x >= 0, m > 1/2, Omega > 0.
%
% EXAMPLE 1:
% % CF of the log-transformed Nakagami RV with m = 2, Omega = 1
%   m     = 2;
%   Omega = 1;
%   t     = linspace(-15,15,501);
%   cf    = cf_LogRV_Nakagami(t,m,Omega);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the log-transformed Nakagami RV with m = 2, Omega = 1')
%
% EXAMPLE 2:
% % PDF/CDF of the log-transformed Nakagami RV with m = 2 Omega = 1
%   m     = 2;
%   Omega = 1;
%   cf    = @(t) cf_LogRV_Nakagami(t,m,Omega);
%   clear options
%   options.N = 2^10;
%   x = linspace(-2,1,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 3: 
% % CF of a linear combination of log-transformed Nakagami RVs
%   m     = [1 1.5 2 3 5];
%   Omega = [1 1   1 3 3];
%   coef  = [1 1   1 1 1];
%   t     = linspace(-5,5,501);
%   cf    = cf_LogRV_Nakagami(t,m,Omega,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of log-transformed Nakagami RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a linear combination of log-transformed Nakagami RVs
%   m     = [1 1.5 2 3 5];
%   Omega = [1 1   1 3 3];
%   coef  = [1 1   1 1 1];
%   cf    = @(t) cf_LogRV_Nakagami(t,m,Omega,coef);
%   clear options
%   options.N = 2^10;
%   x = linspace(-4,4,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% REFERENCES:
% [1] Laurensen, D.I., 1994. Indoor Radio Channel Propagation Modelling by
%     Ray Tracing Techniques.
% [2] Dharmawansa, P., Rajatheva, N. and Ahmed, K., 2007. On the
%     distribution of the sum of Nakagami-m random variables. IEEE 
%     Transactions on Communications, 55(7), 1407-1416.
%
% SEE ALSO:  cf_Nakagami, cf_LogRV_Chi, cf_LogRV_Rice, cf_LogRV_Rayleigh,
% cf_LogRV_MaxwellBoltzmann 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Oct-2018 23:23:19

%% ALGORITHM
%  cf = cf_LogRV_Nakagami(t,m,Omega,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid  = []; end
if nargin < 4, coef  = []; end
if nargin < 3, Omega = []; end
if nargin < 2, m     = []; end

if isempty(m)
    m = 0.5;
end

if isempty(Omega)
    Omega = 1;
end

if isempty(coef) 
    coef = 1;
end

% Check/set equal dimensions for the vectors coef and scale 
[errorcode,coef,Omega,m] = distchck(3,coef(:),Omega(:),m(:));
if errorcode > 0
        error(message('InputSizeMismatch'));
end

% CF of the linear combination of the log-transformed Nakagami RVs
% (expressed by using cf_LogRV_Chi) 
df    = 2*m;
sigma = sqrt(Omega./(2*m));
shift = sum(coef.*log(sigma));
cf    = exp(1i*t*shift) .* cf_LogRV_Chi(t,df,coef,niid);

end
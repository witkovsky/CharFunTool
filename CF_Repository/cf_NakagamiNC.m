function cf = cf_NakagamiNC(t,m,Omega,delta,coef,niid)
%cf_NakagamiNC 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent non-central Nakagami random variables, with the shape
%  parameters m_i > 1/2, spread parameters Omega_i > 0, and the
%  non-centrality parameters delta_i >= 0, for i  = 1,...,N.     
%  
%  The Nakagami distribution (or the Nakagami-m distribution) is a
%  continuous probability distribution for positive-valued random
%  variables, related to the chi-distribution, respectively. 
%
%  In particular, if X ~ NakagamiNC(m,Omega,delta), with the real shape
%  parameter m >= 1/2 and the real spread parameter Omega > 0,  and the
%  non-centrality parameter delta >= 0, Then, the following
%  relationship holds true: X = sqrt(Omega/2*m) * Y, where Y ~
%  Chi(2*m,delta/sqrt(Omega/2*m)).  
%
%  cf_NakagamiNC evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ NakagamiNC(m_i,Omega_i,delta_i)
%  are inedependent non-central Nakagami RVs with the shape parameters m_i
%  >= 1/2, the spread parameters Omega_i > 0, and the non-centrality
%  parameters delta_i >= 0, for i  = 1,...,N.     
%
%  The characteristic function of X ~ NakagamiNC(m,Omega,delta) is defined
%  by 
%   cf_NakagamiNC(t,m,Omega,delta) =
%                  cf_ChiNC((Omega/2*m)*t,df=2*m,delta/sqrt(Omega/2*m)),  
%  where cf_ChiNC(t,df,delta) denotes the characteristic function of the
%  non-central Chi distribution with df degrees of freedom and the
%  non-centrality parameter delta. Hence, the characteristic function of Y
%  is  
%   cf(t) = Prod ( cf_NakagamiNC(t,m_i,Omega_i,delta_i) )
%
% SYNTAX:
%  cf = cf_NakagamiNC(t,m,Omega,delta,coef,niid))
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  m     - vector of the shape parameters of the Nakagami random
%          variables, m > 1/2.  If m is scalar, it is assumed that all
%          shape parameters are equal. If empty, default value is scale =
%          1/2.  
%  Omega - vector of the spread parameters of the Nakagami random
%          variables, Omega > 0.  If Omega is scalar, it is assumed that
%          all spread parameters are equal. If empty, default value is
%          Omega = 1.
%  delta - vector of the non-centrality parameters delta >= 0. If empty,
%          default value is delta = 0. 
%  coef  - vector of the coefficients of the linear combination of the
%          Nakagami random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
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
% % CF of the non-central Nakagami RV with m=2, Omega=1, delta=1
%   m     = 2;
%   Omega = 1;
%   delta = 1;
%   t     = linspace(-10,10,501);
%   cf    = cf_NakagamiNC(t,m,Omega,delta);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the non-central Nakagami RV with m=2, Omega=1, delta=1')
%
% EXAMPLE 2:
% % PDF/CDF of the non-central Nakagami RV with m=2, Omega=1, delta=1
%   m     = 2;
%   Omega = 1;
%   delta = 1;
%   cf    = @(t) cf_NakagamiNC(t,m,Omega,delta);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,3,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 3: 
% % CF of a linear combination of independent non-central Nakagami RVs
%   m     = [1 1 2 3 5];
%   Omega = [1 1 1 3 3];
%   delta = [1 1 1 1 1];
%   coef  = [1 1 1 1 1];
%   t     = linspace(-3,3,501);
%   cf    = cf_NakagamiNC(t,m,Omega,delta,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of non-central Nakagami RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a linear combination of independent non-central Nakagami RVs
%   m     = [1 1 2 3 5];
%   Omega = [1 1 1 3 3];
%   delta = [1 1 1 1 1];
%   coef  = [1 1 1 1 1];
%   cf    = @(t) cf_NakagamiNC(t,m,Omega,delta,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,15,201);
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
% SEE ALSO: cf_Chi

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 06-Oct-2018 00:05:18
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%  cf = cf_NakagamiNC(t,m,Omega,delta,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, niid  = []; end
if nargin < 5, coef  = []; end
if nargin < 4, delta = []; end
if nargin < 3, Omega = []; end
if nargin < 2, m     = []; end

if isempty(m)
    m = 0.5;
end

if isempty(Omega)
    Omega = 1;
end

if isempty(delta)
    delta = 0;
end

if isempty(coef) 
    coef = 1;
end

% Check/set equal dimensions for the vectors coef and scale 
[errorcode,coef,delta,Omega,m] = distchck(4,coef(:),delta(:),Omega(:),m(:));
if errorcode > 0
        error(message('InputSizeMismatch'));
end

% CF of the linear combination of the Nakagami RVs (expressed by using Chi)
df   = 2*m;
coef = coef.*sqrt(Omega./(2*m));
ncp  = delta./sqrt(Omega./(2*m));
cf   = cf_ChiNC(t,df,ncp,coef,niid);

end
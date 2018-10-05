function cf = cf_Chi(t,df,coef,niid)
%cf_Chi 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent CHI-distributed random variables.     
%
%  The CHI-distribution is related to the CHI-SQUARED distribution. It is
%  the distribution of the positive square root of the sum of squares of a
%  set of independent random variables each following a standard normal
%  distribution, or equivalently, the distribution of the Euclidean
%  distance of the random variables from the origin.
%  
%  The most familiar examples are the Rayleigh distribution (which is the
%  CHI-distribution with two degrees of freedom) and the Maxwell�Boltzmann 
%  distribution of the molecular speeds in an ideal gas (CHI-distribution
%  with three degrees of freedom).   
%
%  cf_Chi evaluates the characteristic function cf(t) of Y = sum_{i=1}^N
%  coef_i * X_i, where X_i ~ Chi(df_i) are inedependent CHI-distributed
%  RVs with df_i > 0 degrees of freedom, for i  = 1,...,N.    
%
%  The characteristic function of X ~ Chi(df) is defined by
%   cf_Chi(t,df) = 1F1(df/2,1/2,-t^2/2) + ...
%     + 1i*t*sqrt(2)*gamma((df+1)/2)/gamma(df/2) * 1F1((df+1)/2,3/2,-t^2/2)   
%  Hence, the characteristic function of Y is
%   cf(t) = Prod ( cf_Chi(t,df_i) )
%
% SYNTAX:
%  cf = cf_Chi(t,df,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df    - vector of the degrees of freedom of the chi distributed random
%          variables.  If df is scalar, it is assumed that all degrees of
%          freedom are equal. If empty, default value is df = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          chi distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Chi_distribution
%
% EXAMPLE 1:
% % CF of the distribution of CHI-distributed RV with df = 3
%   df   = 3;
%   coef = 1;
%   t    = linspace(-10,10,501);
%   cf   = cf_Chi(t,df,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('Characteristic function of the \chi-distributed RV with DF=3')
%
% EXAMPLE 2: 
% % CF of a linear combination of independent CHI-distributed RVs
%   df   = [1 1 2 2 3 3 4 4 5 5];
%   coef = [1 1 1 1 1 1 1 1 1 1];
%   t    = linspace(-2,2,501);
%   cf   = cf_Chi(t,df,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent \chi-distributed RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a linear combination of independent CHI-distributed RVs
%   df   = [1 1 2 2 3 3 4 4 5 5];
%   coef = [1 1 1 1 1 1 1 1 1 1];
%   cf   = @(t) cf_Chi(t,df,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(5,25,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 4 (Rayleigh distribution with scale parameter sigma)
% % PDF/CDF of a linear combination of independent Rayleigh RVs
%   df    = [2 2 2 2 2];
%   coef  = [1 1 1 1 1];
%   scale = 2;
%   cf    = @(t) cf_Chi(scale*t,df,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,25,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 5 (Maxwell�Boltzmann distribution with scale parameter sigma)
% % PDF/CDF of a linear combination of independent Maxwell�Boltzmann RVs
%   df    = [3 3 3 3 3];
%   coef  = [1 1 1 1 1];
%   scale = 2;
%   cf    = @(t) cf_Chi(scale*t,df,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(5,30,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% See also: cf_ChiSquare

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 04-Oct-2018 13:47:29

%% ALGORITHM
% cf = cf_Chi(t,df,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 4, niid = []; end
if nargin < 3, coef = []; end
if nargin < 2, df   = []; end

if isempty(df)
    df = 1;
end

if isempty(coef) 
    coef = 1;
end

if isempty(niid)
    niid = 1;
end

% Check/set equal dimensions for the vectors coef and df 
[errorcode,coef,df] = distchck(2,coef(:)',df(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function
szt  = size(t);
t    = t(:);
n    = length(coef);

cf   = 1;
for i = 1:n
    const = sqrt(2) * gamma((df(i)+1)/2) / gamma(df(i)/2);
    cf = cf .* (Hypergeom1F1a(df(i)/2,0.5,-(t*coef(i)).^2/2) + ...
                Hypergeom1F1a((df(i)+1)/2,1.5,-(t*coef(i)).^2/2) .* ...
                (1i*const*t*coef(i)));
end
cf   = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end
function cf = cfX_Gamma(t,alpha,beta)
%cfX_Gamma(t,alpha,beta)  evaluates the characteristic function cf(t) of 
% the GAMMA distribution with the parameters alpha (shape, alpha > 0) and
% beta (rate, beta > 0), i.e. 
%   cf(t) = cfX_Gamma(t,alpha,beta) 
%         = (1 - 1i*t / beta)^(-alpha);
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Gamma_distribution
%
% SYNTAX
%  cf = cfX_Gamma(t,alpha,beta)
%
% EXAMPLE1 (CF of the Gamma distribution with alpha = 2, beta = 2)
%  alpha = 2;
%  beta = 2;
%  t = linspace(-20,20,501);
%  cf = cfX_Gamma(t,alpha,beta);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Gamma distribution with alpha = 2, beta = 2')
%
% EXAMPLE2 (PDF/CDF of the Gamma distribution with alpha = 2, beta = 2)
%  alpha = 2;
%  beta = 2;
%  x = linspace(0,5,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.N = 2^14;
%  cf = @(t) cfX_Gamma(t,alpha,beta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Binomial-Gamma distribution)
%  n = 25;  
%  p = 0.3;
%  alpha = 2;
%  beta = 2;
%  cfX = @(t)  cfX_Gamma(t,alpha,beta);
%  cf = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,25,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
% [1] WITKOVSKY, V.: On the exact computation of the density and of
%     the quantiles of linear combinations of t and F random
%     variables. Journal of Statistical Planning and Inference 94
%     (2001), 1–13.
% [2] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  
% [3] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
%     distribution based on numerical inversion of the compound empirical
%     characteristic function of frequency and severity. Working Paper.
%     Insurance: Mathematics and Economics. 
% [4] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for
%     computing distributions of collective risk models.  Working Paper.
%     Journal of Statistical Software.

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfX_Gamma(t,alpha,beta);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 2, alpha = []; end
if nargin < 3, beta = []; end

%%
if isempty(alpha)
    alpha = 1;
end

if isempty(beta)
    beta = 1;
end

%% Characteristic function of the Gamma distribution
szt = size(t);
t   = t(:);

cf  = (1-1i*t/beta).^(-alpha);
cf  = reshape(cf,szt);
cf(t==0) = 1;

end
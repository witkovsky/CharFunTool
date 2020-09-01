function cf = cfX_PearsonV(t,alpha,beta)
%cfX_PearsonV(t,alpha,beta) Computes the characteristic function of the
% Pearson type V distribution with parameters alpha > 0 (shape) and beta >
% 0 (scale), computed for real vector argument t, i.e. 
%  cf = (2/gamma(alpha)) * (-1i*t/beta).^(alpha/2) .* ...
%        besselk(alpha,2*sqrt(-1i*t/beta));
% where besselk(a,z) denotes the modified Bessel function of the second
% kind. For more details see WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Pearson_distribution
%
% SYNTAX:
%  cf = cfX_PearsonV(t,alpha,beta)
%
% INPUTS:
%  t      -  vector of real values where the CF is evaluated, i.e. CF(t),
%  alpha  -  scalar shape parameter of the Pearson VI distribution, alpha > 0,
%  beta   -  scalar scale parameter of the Pearson VI distribution, beta > 0,
%
% OUTPUT:
% cf     -   calculated values of the characteristic function CF(t)
%
% EXAMPLE1: (CF of the PearsonV distribution)
%  alpha = 3/2;
%  beta = 2/3;
%  t = linspace(-10,10,1001);
%  cf = cfX_PearsonV(t,alpha,beta);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the PearsonV distribution with alpha = 3/2, beta = 2/3')
%
% EXAMPLE2 (PDF/CDF of the PearsonV distribution)
%  alpha = 3/2;
%  beta = 2/3;
%  x = linspace(0,200,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.N = 2^14;
%  options.SixSigmaRule = 10;
%  cf = @(t) cfX_PearsonV(t,alpha,beta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Binomial-PearsonV distribution)
%  n = 25;  
%  p = 0.3;
%  alpha = 3/2;
%  beta = 2/3;
%  cfX = @(t) cfX_PearsonV(t,alpha,beta);
%  cf = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,500,101);
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

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-Sep-2020 13:25:21
%
% Revision history:
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfX_PearsonV(t,alpha,beta);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 3, beta = []; end
if nargin < 2, alpha = []; end
if isempty(beta), beta = 1; end
if isempty(alpha), alpha = 1; end

%% Characteristic function of the PearsonV distribution
szt = size(t);
t   = t(:);

cf = (2/gamma(alpha)) * (-1i*t/beta).^(alpha/2) .* ...
    besselk(alpha,2*sqrt(-1i*t/beta));
cf = reshape(cf,szt);
cf(t==0) = 1;

end

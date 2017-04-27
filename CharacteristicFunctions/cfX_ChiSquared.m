function cf = cfX_ChiSquared(t,df)
%cfX_ChiSquared(t,df) evaluates the characteristic function cf(t) of 
% the CHI-SUQARED distribution with the parameter df (degrees of freedom,
% df >= 1), i.e. 
%   cf(t) = cfX_ChiSquared(t,df) = (1 ? 2?i*?t )^(?df/2) =
%         = cfX_Gamma(t,alpha,beta) = cfX_Gamma(t,df/2,1/2);
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Chi-squared_distribution
%
% SYNTAX
%  cf = cfX_ChiSquared(t,df)
%
% EXAMPLE1 (CF of the ChiSquared distribution with df = 1)
%  df = 1;
%  t = linspace(-50,50,501);
%  cf = cfX_ChiSquared(t,df);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Chi-squared distribution with df = 1')
%
% EXAMPLE2 (PDF/CDF of the ChiSquared distribution with df = 3)
%  df = 3;
%  x = linspace(0,15,101);
%  prob = [0.9 0.95 0.99];
%  cf = @(t) cfX_ChiSquared(t,df);
%  clear options
%  options.xMin = 0;
%  options.xMax = 22;
%  options.N = 2^14;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Binomial-ChiSquared distribution)
%  n = 25;  
%  p = 0.3;
%  df = 3;
%  cfX = @(t) cfX_ChiSquared(t,df);
%  cf = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,80,101);
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
%cf = cfX_ChiSquared(t,df);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);
if nargin < 2, df = []; end

%%
if isempty(df)
    df = 1;
end

%% Characteristic function of the Gamma distribution
szt = size(t);
t   = t(:);

cf  = (1-2i*t).^(-df/2);
cf  = reshape(cf,szt);
cf(t==0) = 1;

end
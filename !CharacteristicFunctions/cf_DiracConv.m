function cf = cf_DiracConv(t,coefs,n_iid)
%%CF_DIRACCONV Characteristic function of the distribution of a linear 
%  combination of independent DIRAC random variables (i.e. nonstochastic
%  constants), Y = sum_{i=1}^N coef_i * Dirac_i, where each Dirac_i
%  represents the constant 1.
%
% SYNTAX
%  cf = cf_DiracConv(t,coefs,n_iid)
%
% EXAMPLE1 (CF of a linear combination of K=100 independent Dirac RVs)
%  K = 100;
%  t = linspace(-10,10,201);
%  idx = 1:K;
%  coefs = 1./((idx - 0.5)*pi).^2;
%  figure; plot(idx,coefs,'.-')
%  title('Coefficients of the linear combination of Dirac RVs')
%  cf = cf_DiracConv(t,coefs);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('Characteristic function of the linear combination of Dirac RVs')
%
% EXAMPLE2 (Compute convolution PDF/CDF of Dirac and Normal RVs by CF2DIST)
%  coef = 5;
%  cf = @(t) cf_DiracConv(t,coef) .* exp(-t.^2/2);
%  result = cf2DIST(cf);

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Mar-2015 10:28:44

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 3, n_iid = []; end
if nargin < 2, coefs = []; end

%%
if isempty(coefs)
    coefs = 1;
end

%% Characteristic function of a linear combination of Dirac variables
smc = sum(coefs);
cf  = exp(1i * t * smc);

if isscalar(coefs) && isscalar(n_iid)
    cf = cf .^ n_iid;
end

end
function cf = cf_chi2Conv(t,coefs,df,ncp,n_iid)
%%CF_CHI2CONV Characteristic function of the distribution of linear 
%  combination of independent chi-squared random variables, i.e. of the
%  distribution of Y = sum_{i=1}^N coef_i * chi2_{df_i,ncp_i}.
%
% SYNTAX
%  cf = cf_chi2Conv(t,coefs,df,ncp,n_iid)
%
% EXAMPLE1 (CF of a linear combination of K=100 independent chi2_1 RVs)
%  N = 2^10;
%  K = 100;
%  t = linspace(-100,100,N);
%  idx = 1:K;
%  coefs = 1./((idx - 0.5)*pi).^2;
%  figure; plot(idx,coefs,'.-')
%  title('Coefficients of the linear combination of \chi^2_1 RVs')
%  df = 1;
%  cf = cf_chi2Conv(t,coefs,df);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('Characteristic function of the linear combination of \chi^2_1 RVs')
%
% EXAMPLE2 (Compute PDF/CDF from the CF by CF2DIST [required])
%  cf = @(t) cf_chi2Conv(t,coefs,df);
%  options.n = N;
%  options.isForcedSymmetric = true;
%  result = cf2DIST(cf,options);
%
% REFERENCES
%  Imhof, J.: Computing the distribution of quadratic forms in normal
%  variables. Biometrika 48, 419–426 (1961).

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Mar-2015 10:28:44

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, n_iid = []; end
if nargin < 4, ncp = []; end
if nargin < 3, df = []; end
if nargin < 2, coefs = []; end

%%
if isempty(ncp) && ~isempty(df)
    isNoncentral = false;
    ncp = 0;
elseif isempty(ncp) && ~isempty(coefs)
    isNoncentral = false;
    ncp = 0;
elseif ~any(ncp)
    isNoncentral = false;
    ncp = 0;
else
     isNoncentral = true;
end

if isempty(df) && ~isempty(coefs)
    df = 1;
elseif isempty(df) && ~isempty(ncp)
    df = 1;
end

if isempty(coefs) && ~isempty(ncp)
    coefs = 1;
elseif isempty(coefs) && ~isempty(df)
    coefs = 1;
end

%% Find the unique coefficients and their multiplicities
if ~isempty(coefs) && isscalar(df) && ~isNoncentral && isempty(n_iid)
    coefs = sort(coefs);
    m    = length(coefs);
    [coefs,idx] = unique(coefs);
    df = df * diff([idx;m+1]);
end

[errorcode,coefs,df,ncp] = distchck(3,coefs(:)',df(:)',ncp(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

% Special treatment for linear combinations with large number of RVs
szcoefs  = size(coefs);
szcoefs  = szcoefs(1)*szcoefs(2);
szt      = size(t);
sz       = szt(1)*szt(2);
szcLimit = ceil(1e3 / (sz/2^16));
idc = 1:fix(szcoefs/szcLimit)+1;

%% Characteristic function of linear combination of noncentral chi-squares
t     = t(:);
idx0  = 1;
cf    = 1;
for j = 1:idc(end)
    idx1 = min(idc(j)*szcLimit,szcoefs);
    idx  = idx0:idx1;
    idx0 = idx1+1;
    aux  = bsxfun(@times,t,coefs(idx));
    if isNoncentral
        aux = bsxfun(@power,(1 -2i * aux),-df(idx)/2) .* ...
            exp(bsxfun(@times,aux,1i*ncp(idx))./(1-2i*aux));
    else
        aux = bsxfun(@power,(1 -2i * aux),-df(idx)/2);
    end
    cf  = cf .* prod(aux,2);
end
cf = reshape(cf,szt);
%cf(t==0) = 1;

if isscalar(coefs) && isscalar(df) && isscalar(ncp) && isscalar(n_iid)
    cf = cf .^ n_iid;
end

end
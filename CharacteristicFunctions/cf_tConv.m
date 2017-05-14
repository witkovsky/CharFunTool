function cf = cf_tConv(t,coefs,df)
%%CF_TCONV Characteristic function of the distribution of a linear
%  combination of independent random variables with Student t-distribution
%  with df degrees of freedom,  i.e. of the distribution of
%  Y = sum_{i=1}^N coef_i * T{df_i}.
%
% SYNTAX
%  cf = cf_tConv(t,coefs,df)
%
% EXAMPLE1 (CF of a linear combination of K=100 independent t RVs)
%  N = 201;
%  K = 50;
%  t = linspace(-10,10,N);
%  idx = 1:K;
%  coefs = 1./((idx - 0.5)*pi).^2;
%  figure; plot(idx,coefs,'.-')
%  title('Coefficients of the linear combination of t RVs')
%  df = 5;
%  cf = cf_tConv(t,coefs,df);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('Characteristic function of the linear combination of t RVs')
%
% EXAMPLE2 (Compute PDF/CDF from the CF by CF2DIST [required])
%  cf = @(t) cf_tConv(t,coefs,df);
%  clear options;
%  options.isForcedSymmetric = false;
%  options.isZeroSymmetric = true;
%  result = cf2DIST(cf,options);
%
% REFERENCES
%  Witkovsky, V.: On the exact computation of the density and of the
%  quantiles of linear combinations of t and F random variables. Journal of
%  Statistical Planning and Inference 94, 1–13 (2001)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 16-Mar-2015 19:38:03

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 3, df = []; end
if nargin < 2, coefs = []; end

%%
if isempty(df) && ~isempty(coefs)
    df = 1;
end

if isempty(coefs) && ~isempty(df)
    coefs = 1;
end

%% Equal size of the parameters
[errorcode,coefs,df] = distchck(2,coefs(:)',df(:)');
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
df2   = df/2;
t     = t(:);
o     = ones(length(t),1);
idx0  = 1;
% cf    = 1;
% for j = 1:idc(end)
%     idx1 = min(idc(j)*szcLimit,szcoefs);
%     idx  = idx0:idx1;
%     idx0 = idx1+1;
%     aux = bsxfun(@times,abs(t),abs(df(idx).*coefs(idx)));
%     aux = bsxfun(@power,aux,df2(idx)) .* ...
%           besselk(ones(length(t),1)*df2(idx),aux);
%     aux = bsxfun(@times,1./(2.^(df2(idx)-1).*gamma(df2(idx))),aux);
%     cf  = cf .* prod(aux,2);
% end

for j = 1:idc(end)
    idx1 = min(idc(j)*szcLimit,szcoefs);
    idx  = idx0:idx1;
    idx0 = idx1+1;
    aux = bsxfun(@times,abs(t),abs(df(idx).*coefs(idx)));
    aux = - aux + bsxfun(@times,df2(idx),log(aux)) + ...
        log(besselk(o*df2(idx),aux,1));
    aux = bsxfun(@plus,aux,(-log(2)*(df2(idx)-1))-gammaln(df2(idx)));
    cf  = exp(sum(aux,2));
end
cf = reshape(cf,szt);
cf(t==0) = 1;

end

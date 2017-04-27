function cf = cf_invGammaConv(t,coefs,A,B,n_iid)
%%CF_INVGAMMACONV Characteristic function of the distribution of a linear 
%  combination of independent inverse Gamma random variables, i.e. of the
%  distribution of Y = sum_{i=1}^N coef_i * IGamma_{A_i,B_i}, where A_i
%  represent the SHAPE parameters and B_i represent the SCALE parameters of
%  the Gamma (!) distribution. For example, IGamma(DF/2,2) represents
%  inverse Chi-squared distribution with DF degrees of freedom.
%
%  Note that there is used also another alternative parametrization, using
%  the reciprocal value(s) of the B parameters (which can be interpreted as
%  a SCALE parameter of the inverse Gamma distribution).
%
% SYNTAX
%  cf = cf_invGammaConv(t,coefs,A,B,n_iid)
%
% EXAMPLE1 (CF of a linear combination of K=100 independent IGamma RVs)
%  N = 2^10;
%  K = 50;
%  t = linspace(-100,100,N);
%  idx = 1:K;
%  coefs = 1./((idx - 0.5)*pi).^2;
%  figure; plot(idx,coefs,'.-')
%  title('Coefficients of the linear combination of IGamma RVs')
%  A = 5/2;
%  B = 2;
%  cf = cf_invGammaConv(t,coefs,A,B);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('Characteristic function of the linear combination of IGamma RVs')
%
% EXAMPLE2 (Compute PDF/CDF from the CF by CF2DIST [required])
%  cf = @(t) cf_invGammaConv(t,coefs,A,B);
%  options.n = N;
%  options.isForcedSymmetric = true;
%  result = cf2DIST(cf,options);
%
% REFERENCES
%  Witkovsky, V.: Computing the distribution of a linear combination of
%  inverted gamma variables. Kybernetika 37, 79–90 (2001).

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Mar-2015 10:28:44

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, n_iid = []; end
if nargin < 4, B = []; end
if nargin < 3, A = []; end
if nargin < 2, coefs = []; end

%%
if isempty(B) && ~isempty(A)
    B = 2;
elseif isempty(B) && ~isempty(coefs)
    B = 2;
elseif ~any(B)
    B = 2;
end

if isempty(A) && ~isempty(coefs)
    A = 2;
elseif isempty(A) && ~isempty(B)
    A = 2;
end

if isempty(coefs) && ~isempty(B)
    coefs = 1;
elseif isempty(coefs) && ~isempty(A)
    coefs = 1;
end

%% Equal size of the parameters   
[errorcode,coefs,A,B] = distchck(3,coefs(:)',A(:)',B(:)');
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
t        = t(:);
idx0     = 1;
cf       = 1;
for j = 1:idc(end)
    idx1 = min(idc(j)*szcLimit,szcoefs);
    idx  = idx0:idx1;
    idx0 = idx1+1;
    aux  = (bsxfun(@times,-1i*t,coefs(idx)./B(idx))).^(1/2);
    aux  = 2*bsxfun(@power,aux,A(idx)) .* ...
           besselk(ones(length(t),1)*A(idx),2*aux);
    aux  = bsxfun(@times,aux,1./gamma(A(idx)));
    cf   = cf .* prod(aux,2);
end
cf = reshape(cf,szt);
cf(t==0) = 1;

if isscalar(coefs) && isscalar(A) && isscalar(B) && isscalar(n_iid)
    cf = cf .^ n_iid;
end

end
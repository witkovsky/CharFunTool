function cf = cf_logBetaConv(t,coefs,A,B)
%%CF_LOGBETACONV Characteristic function of the distribution of linear 
%  combination of independent log-Beta random variables distribution 
%  Y = sum_{i=1}^N coef_i * logBeta_{A_i,B_i}.
%
% SYNTAX
%  cf = cf_logBetaConv(t,coefs,A,B)
%
% EXAMPLE1 (CF of a linear combination of K=100 independent logBeta RVs)
%  N = 2^10;
%  K = 50;
%  t = linspace(-100,100,N);
%  idx = 1:K;
%  coefs = 1./((idx - 0.5)*pi).^2;
%  figure; plot(idx,coefs,'.-')
%  title('Coefficients of the linear combination of IGamma RVs')
%  A = 1;
%  B = 1;
%  cf = cf_logBetaConv(t,coefs,A,B);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('Characteristic function of a linear combination of logBeta RVs')
%
% EXAMPLE2 (Compute PDF/CDF from the CF by CF2DIST [required])
%  cf = @(t) cf_logBetaConv(t,coefs,A,B);
%  options.n = N;
%  options.isForcedSymmetric = true;
%  result = cf2DIST(cf,options);
%
%  REFERENCES:
%  Glaser R.E. (1976). The ratio of the geometric mean and the
%  arithmetic mean for a random sample from a gamma Distribution. Journal
%  of the American Statistical Association, 71(354), 480-487.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Mar-2015 10:28:44

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 4, B = []; end
if nargin < 3, A = []; end
if nargin < 2, coefs = []; end

%%
if isempty(B) && ~isempty(A)
    B = 1;
elseif isempty(B) && ~isempty(coefs)
    B = 1;
end

if isempty(A) && ~isempty(coefs)
    A = 1;
elseif isempty(A) && ~isempty(B)
    A = 1;
end

if isempty(coefs) && ~isempty(B)
    coefs = 1;
elseif isempty(coefs) && ~isempty(A)
    coefs = 1;
end

%% Equal size of the parameters
[errorcode,coefs,A,B] = distchck(3,coefs,A,B);
if errorcode > 0
    error(message('InputSizeMismatch'));
end
A     = A(:);
B     = B(:);
coefs = coefs(:);

%% Characteristic function of linear combination of noncentral chi-squares
aux = 1i*bsxfun(@times,t(:),coefs');
aux = gammalnc(bsxfun(@plus,aux,A')) - ...
      gammalnc(bsxfun(@plus,aux,(A+B)'));
aux = bsxfun(@plus,aux,ones(length(t(:)),1)*(gammalnc(A+B)-gammalnc(A))');
cf  = reshape(prod(exp(aux),2),size(t));
cf(t==0) = 1;

end
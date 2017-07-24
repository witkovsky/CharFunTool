function cf = cf_HAsymptotic(t,alpha,beta,N)
% cf_HAsymptotic Computes the Weibull characteristic function by using the
% the aymptotic expansion (for large values t) of the Fox's H-function for
% 0< beta < 1, see Duby (2017).

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 23-Jul-2017 12:58:35

%% CHECK THE INPUT PARAMETERS
narginchk(3, 4);
if nargin < 4, N = []; end

if isempty(N)
    N = 500;
end

%% ALGORITHM
sz = size(t);
t  = t(:);
z  = -1i*t'*alpha;
n  = (0:N)';

cf = bsxfun(@plus,GammaLog((1+n)*beta)-GammaLog(1+n),-beta*(1+n)*log(z));
cf = sum(bsxfun(@times,beta*(-1).^n,exp(cf)));

cf(z==0) = 1;
cf  = reshape(cf,sz);

end
function cf_conv = cfConv(t,cf,coefs,n)
% CFCONV Characteristic function of a linear combination (convolution) 
%     of iid random variables X_k with given (common) characteristic
%     function CF and the coefficients coefs, such that Y = sum_{k=1}^N
%     coefs(k) * X_k, i.e. cf_conv = Prod_{k=1}^N cf(coefs(k)*t), or
%     cf_conv = cf(t)^N, if all coefs == 1.
%
% SYNTAX:
%     cf_conv = cfConv(t,cf,coefs,n)
%
% INPUT:
%     t         - vector (or array) of input values t where the cf_conv is
%                 evaluated
%     cf        - function handle to the given chracteristic function cf(t)
%     coefs     - vector of coeficients of the linear combination of the 
%                 iid random variables, such that Y = sum(coef(k) * X_k),
%                 and cf_Y = Prod_{k=1}^N cf_X(coef(k)*t).
% 
% OUTPUT:
%     cf_conv   - The characteristic function of a linear combination of
%                 iid RVs with characteristic function cf, evaluated at t.
%
% EXAMPLE:
% % CF of a linear combination of chi-square random variables
%     cf = @(t) cf_AndersonDarlingAsymptotic(t)
%     t = linspace(-10,10,501);
%     coefs = 1./[ 1 2 3 4 5 6 7 8 9 10];
%     cf_conv = @(t) cfConv(t,cf,coefs);
%     plot(t, real(cf_conv(t)),t,imag(cf_conv(t)))
%     clear options;
%     options.nPeriods = 20;
%     options.division1 = [5 8 11 13];
%     options.verbose = false;
%     result = cf2DIST(cf_conv,options)
%
% EXAMPLE: 
% % CF of a linear combination of iid RVs sampled from the empirical
% % distribution function
%     n = 30; p = [0.2 0.7 0.1];
%     data = p(1) * chi2rnd(5,n,1) ...
%            + p(2) * normrnd(10,1,n,1) + p(3) * trnd(1,n,1);
%     clear options;
%     options.SixSigmaRule = 8;
%     options.verbose = false;
%     cf_emp = @(t) cf_DiracMixture(t,data);
%     result_emp = cf2DIST(cf_emp,options)       
%     figure
%     t = linspace(-10,10,501);
%     coefs = 1./[ 1 2 3 4 5 6 7 8 9 10];
%     cf_conv = @(t) cfConv(t,cf_emp,coefs);
%     plot(t, real(cf_conv(t)),t,imag(cf_conv(t)))   
%     result_conv = cf2DIST(cf_conv,options)
%
% DETAILS:
%     A linear combination (convolution) of iid RVs X_k with given
%     characteristic function cf and the coefficients coefs, i.e. that Y =
%     sum_{k=1}^N coefs(k) * X_k, is cf_conv = Prod_{k=1}^N cf(coefs(k)*t),
%     or cf_conv = cf(t)^N, if all coefs == 1.
%
% REFERENCES:
%
% See also cf_t_Conv, cf_chi2_Conv, cf_invGamma_Conv, cf_logBeta_Conv

% Copyright (c) 2015, Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 18-Apr-2015 14:12:32
%
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution
%
% DISCLAIMER THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
% CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
% BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% CHECK THE INPUT PARAMETERS
narginchk(2, 4);
%if nargin < 5, options = []; end
if nargin < 3, coefs = []; end
if nargin < 4, n = []; end

if  isempty(coefs)
    coefs = 1;
end

%% Find the unique coefficients and their multiplicities
if ~isempty(coefs) 
    coefs = sort(coefs);
    m    = length(coefs);
    [coefs,idx] = unique(coefs);    
    df = diff([idx;m+1]);
end

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
idc      = 1:fix(szcoefs/szcLimit)+1;

%% ALGORITHM
t     = t(:);
idx0  = 1;
cf_conv    = 1;
for j = 1:idc(end)
    idx1 = min(idc(j)*szcLimit,szcoefs);
    idx  = idx0:idx1;
    idx0 = idx1+1;
    aux  = bsxfun(@times,t,coefs(idx));
    aux  = bsxfun(@power,cf(aux),df(idx));
    cf_conv  = cf_conv .* prod(aux,2);
end
cf_conv = reshape(cf_conv,szt);
cf_conv(t==0) = 1;

if isscalar(coefs) && isscalar(df) && isscalar(n) 
    cf_conv = cf_conv .^ n;
end

end
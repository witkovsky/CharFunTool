function [qf,result] = cf2QF(cf,prob,options)
%cf2QF  Evaluates the quantiles (QF) from the distribution with 
%      given characteristic function CF by the Gil-Pelaez inversion formula
%      and application of the Newton-Raphson iterative algorithm.  
%
% SYNTAX:
% [qf,result] = cf2QF(cf,prob,options)
% INPUT:
%  cf       - function handle of the characteristic function
%  prob     - vector of probabilities where the QF is computed
%  options  - structure as in cf2CDF_GilPelaez/cf2PDF_GilPelaez
%
% OUTPUT:
%  qf       - vector of the quantile values evaluated at prob.
%  result   - structure with QF further details:
%
% EXAMPLE:
% % QF of the standard Normal distribution
%  cf   = @(t) exp(-t.^2/2);
%  prob = [0.9 0.95 0.99];
%  qf   = cf2QF(cf,prob)
%
% EXAMPLE:
% % QF of the Chi-squared distribution with DF = 1
%  df = 1;
%  cf = @(t) (1 - 2i*t).^(-df/2);
%  prob = [0.9 0.95 0.99];
%  [qf,result] = cf2QF(cf,prob)
%
% EXAMPLE:
% % QF of the linear combination (convolution) of central Chi-squared RVs
%  df   = [1 2 3];
%  cf = @(t) cf_ChiSquare(t,df) ;
%  prob = [0.9 0.95 0.99];
%  [qf,result] = cf2QF(cf,prob)
%
% REFERENCES
%  Imhof, J.: Computing the distribution of quadratic forms in normal
%  variables. Biometrika 48, 419–426 (1961).

% (c) 2015, Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 18-Apr-2015 14:12:32

%% CHECK/SET THE INPUT PARAMETERS
StartTime = cputime;
narginchk(2, 3);
if nargin < 3, options = []; end

if ~isfield(options, 'verbose')
    options.verbose = false;
end

if ~isfield(options, 'qf0')
    options.qf0 = (cf(1e-4)-cf(-1e-4))/(2e-4*1i);
end

if ~isfield(options, 'crit')
    options.crit = 1e-14;
end

if ~isfield(options, 'maxiter')
    options.maxiter = 1000;
end

options.isPlot = false;

%% ALGORITHM
szp  = size(prob);
prob = prob(:);
np   = length(prob);

maxiter   = options.maxiter;
crit      = options.crit;
qf        = options.qf0;
criterion = true;
count     = 0;
while criterion
    count  = count + 1;
    CDFqf  = cf2CDF(cf,qf,options);
    PDFqf  = cf2PDF(cf,qf,options);
    correction  = (CDFqf - prob) ./ PDFqf;
    qf = qf - correction;
    criterion = any(abs(correction) > crit * abs(qf)) ...
        && max(abs(correction)) > crit && count < maxiter;
end

qf    = reshape(qf,szp);
prob  = reshape(prob,szp); 
% Stop the clock
tictoc = cputime - StartTime;

%% RESULTS
result.quantile = qf;
result.prob = prob;
result.CDF = CDFqf;
result.PDF = PDFqf;
result.correction = correction;
result.count = count;
result.cf = cf;
result.options = options;
result.tictoc = tictoc;

end
%%
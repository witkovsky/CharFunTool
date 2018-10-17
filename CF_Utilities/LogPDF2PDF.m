function [x,pdf] = LogPDF2PDF(LogX,LogPDF)
%LogPDF2PDF 
%  Transforms the PDF of the LOG-TRANSFORMED RV (random variable),
%  specified by the known values of LogX and LogPDF, into the PDF of the
%  original RV, specified by the calculated values of x and pdf.
%
% SYNTAX:
%  [x,pdf] = LogPDF2PDF(LogX,LogPDF)
%  [x,pdf] = LogPDF2PDF(result)
%
% INPUT:
%  LogX    - vector of the given x-values of the LOG-TRANSFORMED random
%            variable. Alternatively, it is assumed that it is the given
%            result (structure) from the inversion algorithm (cf2DistGP).
%  LogPDF  - vector of the given PDF values of the LOG-TRANSFORMED random
%            variable. If this parameter is missing,it is assumed that
%            the second parameter is the given result (structure) from the
%            inversion algorithm (cf2DistGP).
%
% OUTPUT:
%  x       - vector of the evaluated x-values of the back transformed
%            random variable, x = exp(LogX).
%  pdf     - vector of the PDF values of the the back transformed
%            random variable, pdf(x) = LogPDF(LogX)/x.
%
% EXAMPLE 1
%  cf = @(t) cf_LogRV_Beta(t,2,3);
%  clear options
%  options.xMax = 0;
%  result  = cf2DistGP(cf,[],[],options);
%  LogX    = result.x;
%  LogPDF  = result.pdf;
%  [x,pdf] = LogPDF2PDF(LogX,LogPDF);
%  figure; plot(x,pdf)
%
% EXAMPLE 2
%  cf = @(t) cf_LogRV_Beta(t,2,3);
%  clear options
%  options.xMax = 0;
%  result  = cf2DistGP(cf,[],[],options);
%  [x,pdf] = LogPDF2PDF(result);
%  figure; plot(x,pdf)
%
% EXAMPLE 3:
%  % PDF of Wilks Lambda (p=5, n=10, q=3) from PDF of LOG-TRANSFORMED RV
%  p    = 5;
%  n    = 10;
%  q    = 3;
%  cf   = @(t) cf_LogRV_WilksLambda(t,p,n,q);
%  clear options
%  options.xMax = 0;
%  result = cf2DistGP(cf,[],[],options);
%  [x,pdf] = LogPDF2PDF(result);
%  figure; plot(x,pdf)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Oct-2018 16:55:59

%% ALGORITHM
%  [x,pdf] = LogPDF2PDF(LogX,LogPDF)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);

if nargin < 2
    if isstruct(LogX)
        result = LogX;
        LogX   = result.x;
        LogPDF = result.pdf;
    else
        error('Missing Inputs')
    end
end

x   = exp(LogX);
pdf = LogPDF./x;

end
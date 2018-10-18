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
%
% EXAMPLE 4:
%  % PDF/CDF/QF/RND of Wilks Lambda (p=5,n=10,q=3) from LOG-TRANSFORMED RV
%  rng(101);
%  p    = 5;
%  n    = 10;
%  q    = 3;
%  cf   = @(t) cf_LogRV_WilksLambda(t,p,n,q);
%  pts  = log(ChebPoints(2^9,[0,1]));
%  clear options
%  options.isPlot = 0;
%  options.xMax = 0;
%  options.xMin = pts(2);
%  xx        = pts(2:end-1);
%  result    = cf2DistGP(cf,xx,[],options);
%  [xxx,cdf] = LogCDF2CDF(result);
%  [~,pdf]   = LogPDF2PDF(result);
%  PDF       = @(x) PDFinterp(x,[0;xxx(:);1],[0;pdf(:);0])
%  CDF       = @(x) CDFinterp(x,[0;xxx(:);1],[0;cdf(:);1])
%  QF        = @(prob) QFinterp(prob,[0;xxx(:);1],[0;cdf(:);1])
%  RND       = @(dim) RNDinterp(dim,[0;xxx(:);1],[0;cdf(:);1])
%  % Plot of PDF/CDF/QF/RND 
%  x   = linspace(0,1,201);
%  p   = linspace(0,1,201);
%  dim = [100000,1];
%  subplot(2,2,1);plot(x,CDF(x))
%  subplot(2,2,2);plot(x,PDF(x))
%  subplot(2,2,3);plot(p,QF(p))
%  subplot(2,2,4);histogram(RND(dim),100,'Normalization','pdf');
%  hold on;plot(x,PDF(x));hold off

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 18-Oct-2018 10:03:08

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
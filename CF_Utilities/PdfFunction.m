function pdf = PdfFunction(xNew,xOld,pdfOld)
%PdfFunction evaluates the interpolant for the PDF at new specified values
%  of x, calculated from the pre-calculated (known) values of xOld, and
%  pdfOld. The evaluation is based on the barycentric interpolation. 
%
% SYNTAX:
%  pdf = PdfFunction(xNew,xOld,pdfOld)
%
% INPUT:
%  xNew    - vector of x-values where the PDF is evaluated, 
%  xOld    - vector of x-values where the values of CDF are known (have
%            been pre-computed). Alternatively, if the this input
%            parameter is missing, it is assumed that the second parameter
%            is the result structure from the cf2DistGP. 
%  cdfOld  - vector of the known CDF values that have been pre-computed 
%            at xOld. If this parameter is missing, it is assumed that the
%            second parameter is the result structure from the cf2DistGP.
%
% EXAMPLE 1
%  cf     = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result = cf2DistGP(cf,[],[],options);
%  xOld   = result.x;
%  pdfOld = result.pdf;
%  pdf = @(x) PdfFunction(x,xOld,pdfOld);
%  x = linspace(-10,10,1001)';
%  plot(x,pdf(x))
%
% EXAMPLE 2
%  cf     = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result = cf2DistGP(cf,[],[],options);
%  pdf = @(x) PdfFunction(x,result);
%  x = linspace(-10,10,1001)';
%  plot(x,pdf(x))

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 16-Oct-2018 13:46:35

%% ALGORITHM
%  pdf = PdfFunction(xNew,xOld,pdfOld)

%% CHECK THE INPUT PARAMETERS
narginchk(2, 3);

if nargin < 3
    if isstruct(xOld)
        result = xOld;
        xOld   = result.x;
        pdfOld = result.pdf;
    else
        error('Missing Inputs')
    end
end

szx  = size(xNew);
xNew = xNew(:);
id   = xNew >= min(xOld) & xNew <= max(xOld);
pdf = zeros(size(xNew));


if any(id)
    pdf(id) = InterpBarycentric(xOld,pdfOld,xNew(id));
end

pdf     = max(0,pdf);
pdf     = reshape(pdf,szx);

end
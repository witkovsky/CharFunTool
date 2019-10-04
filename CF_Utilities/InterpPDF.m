function pdf = InterpPDF(xNew,xGiven,pdfGiven)
%InterpPDF evaluates the interpolant for the PDF at specified values of x,
%  calculated from the pre-calculated (known) values of xGiven, and
%  pdfGiven. The evaluation is based on the barycentric interpolation.
%
% SYNTAX:
%  pdf = InterpPDF(xNew,xGiven,pdfGiven)
%  pdf = InterpPDF(xNew,result)
%
% INPUT:
%  x        - vector of x-values where the PDF is evaluated.
%  xGiven   - vector of x-values where the values of PDF are known (have
%             been pre-computed). Alternatively, if this input parameter is
%             missing, it is assumed that the second parameter is the given
%             result (structure) from the inversion algorithm (cf2DistGP).
%  pdfGiven - vector of the known PDF values that have been pre-computed 
%             at xGiven. If this parameter is missing,it is assumed that
%             the second parameter is the given result (structure) from the
%             inversion algorithm (cf2DistGP).
%
% EXAMPLE 1
%  cf = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result   = cf2DistGP(cf,[],[],options);
%  xGiven   = result.x;
%  pdfGiven = result.pdf;
%  pdf = @(x) InterpPDF(x,xGiven,pdfGiven);
%  x   = linspace(-10,10,1001)';
%  plot(x,pdf(x))
%
% EXAMPLE 2
%  cf = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result = cf2DistGP(cf,[],[],options);
%  pdf    = @(x) InterpPDF(x,result);
%  x      = linspace(-10,10,1001)';
%  plot(x,pdf(x))

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 04-Oct-2019 12:31:36

%% ALGORITHM
%  pdf = InterpPDF(xNew,xGiven,pdfGiven)

%% CHECK THE INPUT PARAMETERS
narginchk(2, 3);

if nargin < 3
    if isstruct(xGiven)
        result = xGiven;
        xGiven   = result.x;
        pdfGiven = result.pdf;
    else
        error('Missing Inputs')
    end
end

szx  = size(xNew);
xNew = xNew(:);
id   = xNew >= min(xGiven) & xNew <= max(xGiven);
pdf = zeros(size(xNew));


if any(id)
    pdf(id) = InterpBarycentric(xGiven,pdfGiven,xNew(id));
end

pdf = max(0,pdf);
pdf = reshape(pdf,szx);

end
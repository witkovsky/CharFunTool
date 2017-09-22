function [pdfI,cdfI] = CFGLSeriesFourierIntegral(t,nCoefs, ...
                       pdfCoefs,cdfCoefs,scale,shift)
% CFGLSeriesFourierIntegral evaluates the Fourier integral of the Legendre
%   series expansion of the PDF and CDF integrand functions (specified by
%   the GAUSS-LEGENDRE SERIES EXPANSION coefficients and the scale and
%   shift parameters from the characteristic function CF) by using the
%   BAKHVALOV-VASILEVA method to evaluate the general Fourier integral,
%      FI(t) = Integral_{A}^B FUN(x) * exp(1i*t*x) dx.
%
% SYNTAX:
%     [pdfI,cdfI] = CFGLSeriesFourierIntegral(t,nCoefs, ...
%                                       pdfCoefs,cdfCoefs,scale,shift)
%
% EXAMPLE: (PDF and CDF by inverting CF of ChiSquared distribution)
%  x   = linspace(1e-12,15,501);
%  cf  = @(t)cf_ChiSquare(t,3);
%  fun = @(t) real(cf(t));
%  A   = 1e-12;
%  B   = 50;
%  [n,pdfCoefs,cdfCoefs,scale,shift] = CFGLSeries(cf,A,B);
%  [I1,I2] = CFGLSeriesFourierIntegral(x,n,pdfCoefs,cdfCoefs,scale,shift);
%  PDF = 2*real(I1)/pi;
%  plot(x,PDF,'.-');grid
%  CDF = 1-2*real(I2)/pi;
%  figure;plot(x,CDF,'.-');grid
%
% REFERENCES:
% [1] BAKHVALOV, N.S., VASILEVA, L.G. Evaluation of the integrals of
%     oscillating functions by interpolation at nodes of Gaussian
%     quadratures. USSR Computational Mathematics and Mathematical Physics,
%     1968, 8(1): 241-249.
% [2] EVANS, G.A.; WEBSTER, J.R. A comparison of some methods for the
%     evaluation of highly oscillatory integrals. Journal of Computational
%     and Applied Mathematics, 1999, 112(1): 55-69.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 21-Sep-2017 09:33:32

%% ALGORITHM CALL
% [pdfI,cdfI] = CFGLSeriesFourierIntegral(t,nCoefs, ...
%                                       pdfCoefs,cdfCoefs,scale,shift)

%% CHECK THE INPUT PARAMETERS
narginchk(1,6);
if nargin < 6, shift    = []; end
if nargin < 5, scale    = []; end
if nargin < 4, cdfCoefs = []; end
if nargin < 3, pdfCoefs = []; end
if nargin < 2, nCoefs   = []; end

if isempty(scale), scale = 1; end
if isempty(shift), shift = 0; end

%% ALGORITHM

szt    = size(t);
t      = t(:)';
id     = t < 0;
t(id)  = -t(id);
K      = scale * exp(1i*t*shift);
t      = scale * t;

if ~isempty(pdfCoefs) && ~isempty(cdfCoefs)
    pdfI     = 0;
    cdfI     = 0;
    for k    = 0:nCoefs-1
        BJ   = 1i^k * (2*k+1) * besselj(k+0.5,t);
        pdfI = pdfI + pdfCoefs(k+1) * BJ;
        cdfI = cdfI + cdfCoefs(k+1) * BJ;
    end
    pdfI     = K .* pdfI ./ sqrt(2*t/pi);
    cdfI     = K .* cdfI ./ sqrt(2*t/pi);
    pdfI(id) = conj(pdfI(id));
    cdfI(id) = conj(cdfI(id));
elseif ~isempty(pdfCoefs) && isempty(cdfCoefs)
    pdfI     = 0;
    for k    = 0:nCoefs-1
        BJ   = 1i^k * (2*k+1) * besselj(k+0.5,t);
        pdfI = pdfI + pdfCoefs(k+1) * BJ;
    end
    pdfI     = K .* pdfI ./ sqrt(2*t/pi);
    pdfI(id) = conj(pdfI(id));
    cdfI     = zeros(size(t));
elseif isempty(pdfCoefs) && ~isempty(cdfCoefs)
    cdfI     = 0;
    for k    = 0:nCoefs-1
        BJ   = 1i^k * (2*k+1) * besselj(k+0.5,t);
        cdfI = cdfI + cdfCoefs(k+1) * BJ;
    end
    cdfI     = K .* cdfI ./ sqrt(2*t/pi);
    cdfI(id) = conj(cdfI(id));
    pdfI     = zeros(size(t));
else
    pdfI     = zeros(size(t));
    cdfI     = pdfI;
end
pdfI         = reshape(pdfI,szt);
cdfI         = reshape(cdfI,szt);

end
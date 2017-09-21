function [I1,I2] = LegendreSeriesFourierIntegral(t,coef1,scale,shift,coef2)
%   LegendreSeriesFourierIntegral evaluates the Fourier integral of the
%   Legendre series expansion of the function FUN (specified by the
%   Legendre coefficients and the scale and shift parameters) by using the
%   BAKHVALOV-VASILEVA method to evaluate  the general Fourier integral,
%     FI(t) = Integral_{A}^B FUN(x) * exp(1i*t*x) dx.
%
% SYNTAX:
%  I = LegendreSeriesFourierIntegral(t,coefs,scale,shift)
%  [I1,I2] = LegendreSeriesFourierIntegral(t,coef1,scale,shift,coef2)
%
% EXAMPLE: (PDF by inverting CF of ChiSquared distribution)
%  x   = linspace(1e-12,15,501);
%  cf  = @(t)cf_ChiSquare(t,3);
%  fun = @(t) real(cf(t));
%  A   = 1e-12;
%  B   = 50;
%  [coefs,scale,shift] = LegendreSeries(fun,A,B);
%  I   = LegendreSeriesFourierIntegral(x,coefs,scale,shift);
%  PDF = 2*real(I)/pi;
%  plot(x,PDF,'.-');grid
%
% EXAMPLE: (PDF by inverting CF of ChiSquared distribution)
%  x       = linspace(1e-12,15,501);
%  cf      = @(t)cf_ChiSquare(t,3);
%  funPdf  = @(t) real(cf(t));
%  funCdf  = @(t) imag(cf(t))./t;
%  A       = 1e-12;
%  B       = 50;
%  [coefs1,scale,shift] = LegendreSeries(funPdf,A,B);
%  coefs2  = LegendreSeries(funCdf,A,B);
%  [I1,I2] = LegendreSeriesFourierIntegral(x,coefs1,scale,shift,coefs2);
%  PDF = 2*real(I1)/pi;
%  CDF = 1-2*real(I2)/pi;
%  plot(x,PDF,'.-');grid
%  figure
%  plot(x,CDF,'.-');grid
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
%[I1,I2] = LegendreSeriesFourierIntegral(t,coef1,scale,shift,coef2)

%% CHECK THE INPUT PARAMETERS
narginchk(1,5);
if nargin < 5, coef2 = []; end
if nargin < 4, shift  = []; end
if nargin < 3, scale  = []; end
if nargin < 2, coef1 = []; end

if isempty(scale), scale = 1; end
if isempty(shift), shift = 0; end

%% ALGORITHM
% Bessel J functions evaluated at required values t
szt    = size(t);
t      = t(:)';
id     = t < 0;
t(id)  = -t(id);
K      = scale * exp(1i*t*shift);
t      = scale * t;

if ~isempty(coef1) && ~isempty(coef2)
    nPts = max(length(coef1)-1,length(coef2)-1);
elseif ~isempty(coef1) && isempty(coef2)
    nPts = length(coef1)-1;
elseif isempty(coef1) && ~isempty(coef2)
    nPts = length(coef2)-1;
else
    nPts = 0;
end

if ~isempty(coef1) && length(coef1) < nPts+1
    coef0 = coef1;
    coef1  = zeros(nPts+1,1);
    coef1(1:length(coef0)) = coef0;
end
if ~isempty(coef2) && length(coef2) < nPts+1
    coef0 = coef2;
    coef2  = zeros(nPts+1,1);
    coef2(1:length(coef0)) = coef0;
end

if ~isempty(coef1) && ~isempty(coef2)
    I1    = 0;
    I2    = 0;
    for k  = 0:nPts
        BJ = 1i^k * (2*k+1) * besselj(k+0.5,t);
        I1 = I1 + coef1(k+1) * BJ;
        I2 = I2 + coef2(k+1) * BJ;
    end
    % Fourier integral evaluated at required values t
    I1     = K .* I1 ./ sqrt(2*t/pi);
    I2     = K .* I2 ./ sqrt(2*t/pi);
    I1(id) = conj(I1(id));
    I2(id) = conj(I2(id));
elseif ~isempty(coef1) && isempty(coef2)
        I1    = 0;
    for k  = 0:nPts
        BJ = 1i^k * (2*k+1) * besselj(k+0.5,t);
        I1 = I1 + coef1(k+1) * BJ;
    end
    I1     = K .* I1 ./ sqrt(2*t/pi);
    I1(id) = conj(I1(id));
    I2     = zeros(size(t));
elseif isempty(coef1) && ~isempty(coef2)
        I2    = 0;
    for k  = 0:nPts
        BJ = 1i^k * (2*k+1) * besselj(k+0.5,t);
        I2 = I2 + coef2(k+1) * BJ;
    end
    I2     = K .* I2 ./ sqrt(2*t/pi);
    I2(id) = conj(I2(id));
    I1     = zeros(size(t));
else         
    I1     = zeros(size(t));
    I2     = I1;
end
I1     = reshape(I1,szt);
I2     = reshape(I2,szt);

end
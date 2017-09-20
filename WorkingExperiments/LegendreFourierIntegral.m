function FI = LegendreFourierIntegral(n,t)
%   LegendreFourierIntegral evaluates the Fourier integrals of the Legendre
%   polynomials used, e.g., in the BAKHVALOV-VASILEVA method to evaluate
%   the general Fourier integral,
%     FI(t) = Integral_{-1}^1 Lk(x) * exp(1i*t*x) dx,
%   for given vector of values omega and all k = 0,...,n.
%
% SYNTAX:
%  FI = LegendreFourierIntegral(n,t)
%
% EXAMPLES:
%  n  = 10;
%  t  = linspace(-5,5,11)';
%  FI = LegendreFourierIntegral(n,t)
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
% Ver.: 20-Sep-2017 22:46:25

%% ALGORITHM
t     = t(:)';
id    = t < 0;
t(id) = -t(id);
nt    = length(t);
K     = scale * exp(1i*t*shift);
t     = scale * t;
FI    = zeros(n+1,nt);
for k = 0:n
    FI(k+1,:) = 1i^k * (2*k+1) * (pi ./ (2*t)).^0.5 .* K .* besselj(k+0.5,t);
end
FI(:,id) = conj(FI(:,id));

end
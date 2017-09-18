function FI = ChebFourierIntegral(n,t)
%   ChebFourierIntegral evaluates the Fourier integrals of the Chebyshev
%   polynomials used, e.g., in the PATTERSON method to evaluate the general
%   Fourier integral, 
%     FI(t) = Integral_{-1}^1 Tk(x) * exp(1i*t*x) dx,
%   for given vector of values omega and all k = 0,...,n.
%
% SYNTAX:
%  FI = ChebFourierIntegral(n,t)
%
% EXAMPLES:
%  n  = 10;
%  t  = linspace(-5,5,11)';
%  FI = ChebFourierIntegral(n,t)
%
% REFERENCES: 
%   EVANS, G.A.; WEBSTER, J.R. A comparison of some methods for the
%   evaluation of highly oscillatory integrals. Journal of Computational
%   and Applied Mathematics, 1999, 112(1): 55-69.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 11-Sep-2017 13:41:45

%% ALGORITHM
t    = t(:)';
nt   = length(t);
idt0 = t == 0;
idt1 = ~idt0;
t    = t(idt1);
FI   = zeros(n+1,nt);
R    = zeros(n+2,n+2);
C    = sqrt(2*pi./t);

if n >= 2
    R(1,1) = 2;
    R(2,2) = 2/3;
    k = 2:2:n; R(3:2:n+1,1) = -2 ./ (k-1) ./ (k+1);
    k = 3:2:n; R(4:2:n+1,2) = -2 ./ (k-2) ./ (k+2);
    for k = 2:2:n
        for j = 0:2:k-1
            R(k+1,j+3) = (k^2-j^2)/(k^2-(j+3)^2)*R(k+1,j+1);
            R(k+2,j+4) = ((k+1)^2-(j+1)^2)/((k+1)^2-((j+1)+3)^2)*R(k+2,j+2);
        end
    end
    for j = 0:n
        FI(j+1,idt0) = R(j+1,1);
        Bj = besselj(j+0.5,t);
        for k = 0:n
            if R(k+1,j+1) ~= 0
                FI(k+1,idt1) = FI(k+1,idt1) + 1i^j*(j+0.5)*R(k+1,j+1) * Bj;
            end
        end
    end
    FI(:,idt1) = bsxfun(@times,FI(:,idt1),C);
end

if n <= 1
    FI(1,idt0) = 2;
    FI(1,idt1) = C .* besselj(0.5,t);
end

if n == 1
    FI(2,idt0) = 0;
    FI(2,idt1) = 1i * C .* besselj(1.5,t);
end   
    
end
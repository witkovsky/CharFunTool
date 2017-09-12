function I = ChebIntegral(omega,nPts)
%   ChebIntegral evaluates the required Fourier integrals of Chebyshev
%   polynomials used in the Pattersom Fourier integral method,
%     I(omega) = Integral_{-1}^1 Tk(x) * exp(1i*omega*x) dx,
%   for given vector of values omega and all k = 0,...,nPts.
%
% SYNTAX:
%  I = ChebIntegral(omega,nPts)
%
% EXAMPLES:
%  omega = linspace(-5,5)';
%  nPts  = 10;
%  I = ChebIntegral(omega,nPts)
%
% REFERENCES: 
%   EVANS, G.A.; WEBSTER, J.R. A comparison of some methods for the
%   evaluation of highly oscillatory integrals. Journal of Computational
%   and Applied Mathematics, 1999, 112(1): 55-69.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 11-Sep-2017 13:41:45

%% ALGORITHM

omega   = omega(:)';
nOmega  = length(omega);
I       = zeros(nPts+1,nOmega);
kE      = 0:2:nPts;
kO      = 1:2:nPts;
R       = zeros(1,nPts+1);
R(kE+1) = ((-1).^kE+1)./(1-kE.^2);
R(kO+1) = -2./((kO-2).*(kO+2));
R(2)    = 2/3;

for k  = 1:2:(nPts+1)
    J  = 0;
    Rk = R(k); 
    for j  = 1:2:k
        Rk = Rk *((k-1)^2-(j-1)^2) / ((k-1)^2-((j-1)+3)^2);       
        c  = Rk * (2*(j-1)+1) / 2;
        J  = J + 1i^(j-1) * c * besselj((j-1)+0.5,omega);
    end
    I(k,:) = sqrt(2*pi./omega) .* J;    
end

for k  = 2:2:(nPts+1)
    J  = 0;
    Rk = R(k);
    for j  = 2:2:k
        Rk = Rk*((k-1)^2 - (j-1)^2)/((k-1)^2-((j-1)+3)^2);      
        c  = Rk * (2*(j-1)+1)/2;
        J  = J + 1i^(j-1) * c * besselj((j-1)+0.5,omega);
    end
    I(k,:) = sqrt(2*pi./omega) .* J;     
end
I(:,omega==0) = 0;

end
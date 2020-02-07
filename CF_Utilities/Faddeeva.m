function w = Faddeeva(z)
% Faddeeva 
%  Computes the Faddeeva function or the Kramp function, which is a scaled
%  complex complementary error function, in complex argument z. The
%  algorithm Faddeeva was adapted from the original algorithm fadf.m
%  written by Sanjar M. Abrarov and Brendan M. Quine.
%
%  The Faddeeva function is defined as
%   w(z) = exp(-z^2) * erfc(-1i*z) 
%        = exp(-z^2) * (1 - erf(-1i*z))
%        = exp(-z^2) * (1 + 2i/sqrt(pi) * int_0^z exp(t^2)dt),
%  where erf(z) is the complex error function and erfc(z) is the complex
%  complementary error function.
%
%  In particular, the following relations hold true:
%   erfcx(z)  = exp(z^2) * erfc(z)
%             = w(1i*z)
%   erfc(z)   = exp(-z^2) * w(1i*z) 
%   erf(z)    = 1 - erfc(z)
%             = 1 - exp(-z^2) * w(1i*z)
%   erfi(z)   = -1i * erf(1i*z)
%             = -1i * (1 - exp(-z^2) * w(1i*z))
%   Dawson(z) = sqrt(pi)/2 * exp(-z^2) * erfi(z)
%             = 1i * sqrt(pi)/2 * (exp(-z^2) - w(z)).
%
%  The algorithmic implementation utilizes the approximations based on the
%  Fourier expansion and the Laplace continued fraction. The code covers
%  with high-accuracy the entire complex plain required for practical
%  applications. The code remains highly accurate even at vanishing
%  imaginary argument y -> 0, where y = Im[z]. The worst detected accuracy
%  is ~10^-13.
%
%  For more details see also
%  http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
%
%  SYNTAX:
%   w = Faddeeva(z)
%
%  INPUTS:
%   z      - possibly complex input argument (vector or array),
%
%  OUTPUTS:
%   w      - values of the Faddeeva function
%
%  EXAMPLE 1:
%  % Faddeeva function
%  z = 1i + linspace(-5,5);
%  w = Faddeeva(z);
%  plot(w)
%
%  EXAMPLE 2:
%  % erfc(z) is the complex complementary error function.
%  z    = 1i + linspace(-5,5);
%  erfc = @(z) exp(-z.^2) .* Faddeeva(1i*z);
%  f    = erfc(z);
%  plot(f)
%
%  EXAMPLE 3:
%  % erf(z) is the complex error function.
%  z    = 1i + linspace(-5,5);
%  erf  = @(z) 1 - exp(-z.^2) .* Faddeeva(1i*z);
%  f    = erf(z);
%  plot(f)
%
%  EXAMPLE 4:
%  % erfi(z) is the imaginari complex error function.
%  z    = 1i + linspace(-5,5);
%  erfi  = @(z) -1i * (1 - exp(z.^2) .* Faddeeva(-z));
%  f    = erfi(z);
%  plot(f)
%
%  EXAMPLE 5:
%  % Dawson(z) is the Dawson complex function.
%  z    = 1i + linspace(-5,5);
%  Dawson  = @(z) 1i * sqrt(pi)/2 * (exp(-z.^2) - Faddeeva(z)));
%  f    = Dawson(z);
%  plot(f)
%
%  REFERENCES:
%  [1] S. M. Abrarov and B. M. Quine, Appl. Math. Comput., Efficient
%      algorithmic implementation of the Voigt/complex error function based
%      on exponential series approximation, 218 (2011) 1894-1902.
%      http://doi.org/10.1016/j.amc.2011.06.072
%  [2] S. M. Abrarov and B. M. Quine, On the Fourier expansion method
%      for highly accurate computation of the Voigt/complex error function
%      in a rapid algorithm, arXiv:1205.1768v1 (2012).
%      http://arxiv.org/abs/1205.1768
%  [3] W. Gautschi, Efficient computation of the complex error function,
%      SIAM J. Numer. Anal., 7 (1970) 187-198.
%      http://www.jstor.org/stable/2949591
%  [4] T. M. Karbach, G. Raven and M. Schiller, Decay time integrals in
%      neutral meson mixing and their efficient evaluation,
%      arXiv:1407.0748v1 (2014).
%      http://arxiv.org/abs/1407.0748
%
%  AUTHOR/CREDITS:
%   This is adapted version of the original algorithm fadf.m  written by
%   Sanjar M. Abrarov and Brendan M. Quine, York University, Canada. For
%   more details see the LICENCE.
%   github.com/aoeftiger/faddeevas/blob/master/cernlib_root_extended/fadf.m
%
%  LICENCE:
%  Copyright (c) 2016, Sanjar Abrarov
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are
%  met:
%
%      * Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in
%        the documentation and/or other materials provided with the
%        distribution 
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
%  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
%  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
%  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
%  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
%  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
%  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
%  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This is modified version of the original MATLAB code fadf.m by Sanjar M.
% Abrarov and Brendan M. Quine
% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Sep-2019 20:26:55

%% FUNCTION CALL
% w = Faddeeva(z)

%% CHECK THE INPUT PARAMETERS
narginchk(1,1);
sz = size(z);
z  = z(:);

%% ALGORITHM
ind_neg    = imag(z)<0;        % if some imag(z) values are negative, then
z(ind_neg) = conj(z(ind_neg)); % bring them to the upper-half plane
ind_ext    = abs(z)>15;        % external indices
ind_band   = abs(z)<=15 & imag(z)<10^-5; % narrow band indices
w          = zeros(size(z));   % define array

% internal area. This area is the most difficult for accurate computation
w(~ind_ext & ~ind_band) = fexp(z(~ind_ext & ~ind_band));
w(ind_ext)              = contfr(z(ind_ext)); % external area
w(ind_band)             = smallim(z(ind_band)); % narrow band

% Convert for negative imag(z) values
w(ind_neg) = conj(2*exp(-z(ind_neg).^2) - w(ind_neg));
w = reshape(w,sz);
end
%% FUNCTION FEXP
function FE = fexp(z,tauM) 
% FEXP - the Fourier expansion approximation

if nargin == 1
    tauM = 12; % default margin value
end
maxN = 23; % number of summation terms

n  = 1:maxN;
aN = 2*sqrt(pi)/tauM*exp(-n.^2*pi^2/tauM^2); % Fourier coefficients

z1 = exp(1i*tauM*z); % define first repeating array
z2 = tauM^2*z.^2;    % define second repeating array

FE = sqrt(pi)/tauM*(1 - z1)./z2; % initiate array FE
for n = 1:maxN
    FE = FE + (aN(n)*((-1)^n*z1 - 1)./(n^2*pi^2 - z2));
end
FE = 1i*tauM^2*z/sqrt(pi).*FE;
end

%% FUNCTION CONTFR
function CF = contfr(z) 
% CONTFR - the Laplace continued fraction approximation

aN = 8; % initial integer
aN = 1:2*aN;
aN = aN/2;

CF = aN(end)./z; % start computing from the last aN
for n = 1:length(aN) - 1
    CF = aN(end-n)./(z - CF);
end
CF = 1i/sqrt(pi)./(z - CF);
end

%% FUNCTION SMALLIM
function SIm = smallim(z) 
% SMALLIM - approximation at small imag(z)

ind_0 = abs(real(z))<=1e-3; % indices near origin
ind_4 = abs(real(z))>4; % indices for |Re[z]| > 4

% If 10^-3 < |Re[z]| <= 4, then:
SIm(~ind_0 & ~ind_4) = fexp(z(~ind_0 & ~ind_4),12.1);

% else if |Re[z]| <= 10^-3, then:
SIm(ind_0) = (1-z(ind_0).^2).*(1+2i/sqrt(pi)*z(ind_0));

del    = 10^-5; % assign delta
z_del  = real(z(ind_4)) + 1i*del; % incorporate delta
zr     = imag(z(ind_4))./del; % assign ratio
ff_ref = fexp(z_del); % assign reference of the Faddeeva function

% Finally, for |Re[z]| > 4 apply
SIm(ind_4) = zr.*real(ff_ref) + (1 - zr).*exp(-real(z_del).^2) ...
    + 1i*imag(ff_ref);
end
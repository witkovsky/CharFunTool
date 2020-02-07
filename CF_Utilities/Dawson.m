function D = Dawson(z)
%Dawson 
%  Computes the Dawson function for complex input z (scalar, vector, or
%  array). Here, the Dawson function is defined by 
%   D(z) = exp(-z^2) * integral_0^{z} exp(t^2) dt
%        = integral_0^1 exp(z^2*(z^2-1)^2) dz
%        = 1/2 * integral_0^Inf exp(-t^2/2)*sin(z*t) dt
%  The Dawson function is the one-sided Fourier–Laplace sine transform of
%  the Gaussian function. It is closely related to the error function erf:
%   D(z) = sqrt(pi)/2 * exp(-z^2) * erfi(z),
%        = 1i * sqrt(pi)/2 * (exp(-z^2) - Faddeeva(z)).
%  where erfi(z) is the imaginary error function, erfi(z) = 1i*erf(1i*z)
%  and Faddeeva(z) is the Faddeeva function or the Kramp function, which is
%  a scaled complex complementary error function, in complex argument z.
%  For real x, 
%   D(x) = sqrt(pi)/2 * imag(Faddeeva(x)).
%
% SYNTAX
%  D = Dawson(z)
%
% INPUTS:
%  z     - vector or array of (complex) values, where the Dawson function 
%          is evaluated. 
%
% OUTPUT:
%  D     - vector or array of the calculated Dawson function values.
%  
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Dawson_function
%
% EXAMPLE 1:
%   x = linspace(-20,20,501);
%   D = Dawson(x);
%   plot(x,D);grid
%   xlabel('x')
%   title('Dawson function')
%
% EXAMPLE 2 (CF of the Rayleigh distribution)
%  % CF of the Rayleigh distribution with the scale parameter sigma = 1
%   scale = 1;
%   cf    = @(t) (1-2*(t*scale/sqrt(2)) .* Dawson((t*scale/sqrt(2)))) + ...
%           (1i*scale*t*sqrt(2)*gamma(3/2)) .* exp(-(t*scale/sqrt(2)).^2);
%   t     = linspace(-20,20,501);
%   plot(t,real(cf(t)),t,imag(cf(t)));grid
%   xlabel('t')
%   ylabel('CF')
%   title('CF of the Rayleigh distribution with the parameter sigma = 1')
%
% EXAMPLE 3:
%  % CDF/PDF of the Rayleigh distribution with the scale parameter sigma = 1
%   scale = 1;
%   cf    = @(t) (1-2*(t*scale/sqrt(2)) .* Dawson((t*scale/sqrt(2)))) + ...
%           (1i*scale*t*sqrt(2)*gamma(3/2)) .* exp(-(t*scale/sqrt(2)).^2);
%   x = linspace(0,5);
%   prob = [0.9 0.95 0.99];
%   clear options;
%   options.N = 2^10;
%   result = cf2DistGP(cf,x,prob,options)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Sep-2019 20:26:55

%% ALGORITHM

sz = size(z);
z  = z(:);

if isreal(z)
    D = (sqrt(pi)/2) * imag(Faddeeva(z));
else
    D = (1i * sqrt(pi)/2) * (exp(-z.^2) - Faddeeva(z));
end

D = reshape(D,sz);
end
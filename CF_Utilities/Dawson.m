function F = Dawson(x,tol)
%Dawson 
%  Computes the Dawson function for complex input x (scalar, vector, or
%  array). Here, the Dawson function is defined by 
%   F(x) = exp(-x^2) * integral_0^{x} exp(y^2) dy
%        = integral_0^1 exp(x^2*(z^2-1)^2) dz
%        = 1/2 * integral_0^Infty exp(-t^2/2)*sin(x*t) dt
%  The Dawson function is the one-sided Fourier–Laplace sine transform of
%  the Gaussian function. It is closely related to the error function erf:
%   F(x) = sqrt(pi)/2 * exp(-x^2)  * erfi(x),
%  where erfi is the imaginary error function, erfi(x) = 1i*erf(1i*x). 
%
% SYNTAX
%  F = Dawson(x,tol)
%
% INPUTS:
%  x     - vector or array of (complex) values, where the Dowson function F
%          is evaluated. 
%  tol   - parameter of the relative tolerance RelTol used for integration.
%          If empty, default value is tol = 1e-6.  
%
% OUTPUT:
%  F     - vector or array of the calculated function values.
%  
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Dawson_function
%
% EXAMPLE 1:
%   x = linspace(-20,20,501);
%   F = Dawson(x);
%   plot(x,F);grid
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
% Ver.: 09-Oct-2018 14:54:14

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);

if nargin < 2
    tol = []; 
end

if isempty(tol)
    tol = 1e-6;
end

sz = size(x);
x  = x(:);
F  = x .* integral(@(z) exp((z.^2-1).*x.^2),0,1,'ArrayValued',true, ...
    'RelTol',tol);
F  = reshape(F,sz);

end
function p = PhiZX(z)
%% PhiZX
%  Computes the CDF for the standard normal distribution in complex
%  argument z by using the Faddeeva function.
%
% CDF for the normal distribution.
%  EXAMPLE 1:
%   z = 1i + linspace(-5,5)';
%   f = Phi(z);
%   disp([z f])

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Jan-2025 12:04:15

p = 0.5 * erfcZX(-z / sqrt(2));

end
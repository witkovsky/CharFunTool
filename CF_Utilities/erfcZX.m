function f = erfcZX(z)
% erfcZX
%  Computes the complementary error function in complex argument z by using
%  the Faddeeva function.
%
%  The erfcZX function is defined as
%   erfcZX(z) = exp(-z^2) * w(1i*z)
%  where w(z) is the Faddeeva function or the Kramp function, which is a
%  scaled complex complementary error function, in complex argument z.
%
%  SYNTAX:
%   f = erfcZX(z)
%
%  INPUTS:
%   z      - possibly complex input argument (vector or array),
%
%  OUTPUTS:
%   f      - values of the erfcZX function
%
%  EXAMPLE 1:
%  % erfcZX the complementary error function in complex argument z
%  z = 1i + linspace(-5,5);
%  f = erfcZX(z);
%  plot(f)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Sep-2019 20:26:55
%% FUNCTION CALL
% f = erfcZX(z)

%% ALGORITHM
if isreal(z)
    f = erfc(z);
else
    f = exp(-z.^2) .* Faddeeva(1i*z);
end
end
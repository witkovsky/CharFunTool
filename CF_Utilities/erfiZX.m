function f = erfiZX(z)
% erfiZX
%  Computes the imaginary error function in complex argument z by using the
%  Faddeeva function.
%
%  The erfiZX function is defined as
%   erfiZX(z) = -1i * erfZX(1i*z)
%             = -1i * (1 - exp(z^2) * w(-z))
%  where w(z) is the Faddeeva function or the Kramp function, which is a
%  scaled complex complementary error function, in complex argument z.
%
%  SYNTAX:
%   f = erfiZX(z)
%
%  INPUTS:
%   z      - possibly complex input argument (vector or array),
%
%  OUTPUTS:
%   f      - values of the erfiZX function
%
%  EXAMPLE 1:
%  % erfiZX the imaginary error function in complex argument z
%  z = 1i + linspace(-5,5);
%  f = erfiZX(z);
%  plot(f)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Sep-2019 20:26:55
%% FUNCTION CALL
% f = erfiZX(z)

%% ALGORITHM

f = -1i * erfZX(1i*z);

end
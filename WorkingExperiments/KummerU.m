function [y,details] = KummerU(a,b,z)
%KummerU Computes U(a,b,z) - the Kummer's confluent hypergeometric function 
% of the second kind, for complex arguments a,b,z, where Re(a)>0, Re(z)>=0.
%
% KummerU is based on integral representation of U(a,b,z), which is
% evaluated by numerical integration, U(a,b,z) = int_0^inf fun(a,b,z,t) dt,
% for properly chosen integrand function fun(a,b,z,t).
%
% In this implementation, the integral is numerically evaluated by the
% MATLAB function INTEGRAL, i.e. Q = INTEGRAL(FUN,A,B) approximates the
% integral of function FUN from A to B using global adaptive quadrature and
% default error tolerances.
%
% The standard integral representation is formally correct for Re(a)>0,
% Re(z)>0. However, KummerU computes correct values also for z some values
% from the region with Re(z)<=0. 
% In particular, KummerU gives correct results for the following situations
% 1. For all z such that Re(z) = 0 (PURLY IMAGINARY z), and Re(a) > 0.
% 2. For Re(z) < 0 & abs(Im(z)/Re(z)) > 100 (LARGE OSCILATIONS), Re(a) > 0.
% 3  For Re(z) < 0 & abs(Im(z)/Re(z)) < 100 (SMALL OSCILATIONS), Re(a) > 0
%    and  Re(a-b+1)>0.
%
% TODO (The situations which are not covered by KummerU):
% 1. Implementation for Re(a) <= 0 and arbitrary b, z;
% 2. Implementation for Re(a) > 0 & Re(a-b+1)<0 and Re(z) < 0 with 
%    abs(Im(z)/Re(z)) < 100.
%
% SYNTAX
%   [y,details] = KummerU(a,b,z)
%
% EXAMPLE 1
% a = 3;
% b = 2;
% z = 1 + 1i*(0:0.1:1)';
% [y,details] = KummerU(a,b,z)
%
% EXAMPLE 2 (CF of the Fisher-Snedecor F(d1,d2) distribution)
% % See e.g. https://en.wikipedia.org/wiki/F-distribution
% d1 = 5;
% d2 = 2;
% t = linspace(-20,20,501)';
% cf = gamma((d1+d2)/2) *  KummerU(d1/2,1-d2/2,-(d2/d1)*1i*t) / gamma(d2/2);
% figure; plot(t,real(cf),t,imag(cf))
% title('CF of the Fisher-Snedecor F(5,2) distribution')
% xlabel('t')
% ylabel('CF')

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Mar-2016 12:34:01

%% CHECK THE INPUT PARAMETERS
narginchk(2, 3);
if nargin < 3, z = []; end

if isempty(z)
    y = [];
    details.title   = 'Integrand FUN of Confluent Hypergeometric Function: y = U(a,b,z)';
    details.fun     = @(a,b,z,t)fun(a,b,z,t);
    return
end

if real(a) <= 0
    w = 'Integral representation of U(a,b,z) requires Re(a) > 0';
    error(w);
end

%% EVALUATE the function U(a,b,z)
sz = size(z);
z = z(:);
n = length(z);
y = zeros(n,1);
for k=1:n
%     zr = real(z(k));
%     zi = imag(z(k));
%     if (zi ~= 0 && zr/abs(zi) < -100)
%         w = 'Integral representation of U(a,b,z) requires Re(z) > 0';
%         warning(w)
%     end
    y(k) = integral(@(t)fun(a,b,z(k),t),0,inf);
end
y = reshape(y,sz);

%% SET THE RESULTS
if nargout > 1
    details.title   = 'Confluent Hypergeometric Function: y = U(a,b,z)';
    details.y       = y;
    details.a       = a;
    details.b       = b;
    details.z       = z;
    details.fun     = @(z,t)fun(a,b,z,t);
end

end
%% Function integrand
function f = fun(a,b,z,t)
%FUN Evaluates the integrand of the Kummer's Confluent Hypergeometric
%    function of the second kind U(a,b,z), for Re(a) > 0. Here,
%    U(a,b,z) = Int_0^inf fun(a,b,z,t) dt

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Mar-2016 12:34:01

zr = real(z);
zi = imag(z);
% PROBLEM: This assumes that Re(a) > 0
% TODO: Find alternative solution for Re(a) <=0
if (zr == 0 && zi == 0) % checked OK /WELL DEFINED
    f = t.^(a-1) .* (1+t).^(b-a-1) / gamma(a);
elseif (zr == 0 && zi ~= 0) % checked OK /WELL DEFINED
    if (abs(zi) >= 1)
        % 1st order FIT / for LARGE OSCILLATIONS
        f = (-1i/(zi*gamma(a))) .*  (-(1i/zi)*t).^(a-1) .* ...
            (1-(1i/zi)*t).^(b-a-1) .* exp(-t);
    elseif (abs(zi) < 1)  % 1st order FIT / for SMALL OSCILLATIONS       
        if (zi > 0)
        % Change of variables for stabilizing the integration domain (0,T)
        f = (-1i/gamma(a)) .*  (-1i*t).^(a-1) .* (1-1i*t).^(b-a-1) .* exp(-zi*t);
        elseif (zi < 0)
            % Change of variables for stabilizing the integration domain (0,T)
            f = (1i/gamma(a)) .*  (1i*t).^(a-1) .* (1+1i*t).^(b-a-1) .* exp(zi*t);
        end
    end
elseif (zr ~= 0 && zi == 0) % checked OK /WELL DEFINED
    if (abs(zr) > 100 || zr > 0)
        % OK /WELL DEFINED (abs(zr) > 100 || zr > 0),
        % f = exp(-zr*t) .* t.^(a-1) .* (1+t).^(b-a-1) / gamma(a);
        % Change of variables for stabilizing the integration domain (0,T)
        % WELL DEFINED for zr > 0, but also for large negative zr,
        % such that 1+T/zr ~= 0 and T < -zr, say zr for zr < -100
        f = (1/zr) * exp(-t) .* (t/zr).^(a-1) .* (1+t/zr).^(b-a-1) / gamma(a);
    elseif (zr >= -100) 
        % assumes Re(a-b+1)>0,
        % OK for SMALL negative z, say -100 < z < 0,
        f = (2/gamma(a)/gamma(a-b+1)) * z^(1/2-b/2) * exp(-t) .* ...
            t.^(a-b/2-1/2) .* besselk(b-1,2*sqrt(z*t));
    end
elseif (zr ~= 0 && zi ~= 0)  % checked OK /WELL DEFINED
    if (abs(zr/zi) > 1) % = abs(zi/zr) < 1 => SMALL OSCILLATIONS
        if (zr > 0)
            % 2nd order FIT (Fourier Integral Transform):
            % OK / WELL DEFINED for SMALL OSCILLATIONS, zr > 0 & abs(zi/zr) < 1
            f = (1/(zr*gamma(a))) .* (t/zr).^(a-1) .* (1+(t/zr)).^(b-a-1) .* ...
                exp(-(1+1i/(zr/zi))*t);
            % THE ORIGINAL INTEGRAL FORM is also WELL DEFINED change of
            % variables stabilize the effective integration domain
            %f = (1/gamma(a)) .* (t).^(a-1) .* (1+t).^(b-a-1) .* exp(-z*t);
        elseif (zr < 0)
            if (abs(zr) > 100)
                % f = (-1/gamma(a)) .* (-t).^(a-1) .* (1-t).^(b-a-1) .* exp(z*t);
                f = (1/(zr*gamma(a))) .* (t/zr).^(a-1) .* (1+t/zr).^(b-a-1) .* ...
                    exp(-(1+1i/(zr/zi))*t);
            elseif (abs(zr) <= 100)
                f = (2/gamma(a)/gamma(a-b+1)) * z^(1/2-b/2) * exp(-t) .* ...
                    t.^(a-b/2-1/2) .* besselk(b-1,2*sqrt(z*t));
            end
        end
    elseif (abs(zr/zi) <= 1) % (zr ~= 0 & abs(zi/zr) > 1), i.e. LARGE OSCILLATIONS
        % 1st order FIT / suggested for LARGE OSCILLATIONS, abs(zi/zr) > 1
        % OK / WELL DEFINED for SMALL OSCILLATIONS, zr > 0 & abs(zi/zr) < 1
        % PROBLEM: Oscillatory fun for large abs(zr/zi)> 1
        % It works OK also for NEGATIVE zr with LARGE OSCILLATIONS
        % (zr < 0 & abs(zi/zr) > 1)
        f = (-1i/(zi*gamma(a))) .* exp(1i*(zr/zi)*t) .* ...
            (-(1i/zi)*t).^(a-1) .* (1-(1i/zi)*t).^(b-a-1) .* exp(-t);
    end
end

end

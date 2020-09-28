function [B,D,A,c] = PolyCalUpdate(mu0,nu0,beta0,m,order,q,Aijk)
% PolyCalUpdate is an auxiliary function for the algorithm PolyCal which
% calculates the update of the matrices A, B and the vector c,
% for a polynomial comparative calibration model of given order, locally at
% the specified mu0, nu0, and beta0. If the coefficient
% (m,order+1,q)-dimensional array Aijk is given, then the PolyCalUpdated
% matrices A, B and the vector c are calculated for the generalized
% polynomial calibration model.
%
% SYNTAX
%  [B,D,A,c] = PolyCalUpdate(mu0,nu0,beta0,n,order,q,Aijk)
%
% INPUTS:
%  mu0   - m-dimensional vector of initial values of the vector parameter
%          mu (mean value of the measurement vector X) used for
%          linrarization of the calibration model.
%  nu0   - m-dimensional vector of initial values of the vector parameter
%          nu (mean of the measurement vector Y) used for linrarization of
%          the calibration model.
%  beta0 - (order+1)-dimensional vector of initial values of the vector
%          parameter beta (parameters of the polynomial calibration
%          function of given order) used for linrarization of the
%          calibration model.
%  m     - dimension of the parameter vectors mu0 and nu0.
%  order - order of the polynomial calibration function. Hence, the
%          dimension of the parameter vector beta0 is (order+1).
%  q     - third dimension of the (m,order+1,q)-dimensional array of
%          coefficients Aijk, which specifies the generalized polynomial
%          calibration model.
%  Aijk  - the (m,order+1,q)-dimensional array of coefficients Aijk, which
%          specifies the generalized polynomial calibration model.
%
% EXAMPLE
%  mu0   = [-19.6843   -9.8142    0.0989   10.0149   19.9634   29.8131]';
%  nu0   = [92.1489   96.0965  100.0499  103.9924  107.9354  111.8316]';
%  beta0    = [100.0108    0.3982   -0.0001]';
%  m     = 6;
%  order = 2;
%  q     = [];
%  Aijk  = [];
% [B,D,A,c] = PolyCalUpdate(mu0,nu0,beta0,m,order,q,Aijk)
%
% REFERENCES:
% [1] WITKOVSKY V. and WIMMER G.: PolyCal - MATLAB algorithm for
% comparative polynomial calibration and its applications. AMCTM 2020.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-Sep-2020 16:46:23

%% Algorithm
if nargin < 7
    Aijk = [];
end

if nargin < 6
    q = [];
end

polyPars = flip(beta0);
if isempty(Aijk)
    polyParsD1 = polyder(polyPars);
    D = diag(polyval(polyParsD1,mu0));
    A = [D -eye(m)];
    B = ones(m,order+1);
    for j = 1:order
        B(:,j+1) = mu0.*B(:,j);
    end
    c = polyval(polyPars,mu0) - nu0;
else
    B = zeros(m,order+1);
    D = zeros(m,order+1);
    for i = 1:m
        for j = 1:(order+1)
            for k = 1:q
                B(i,j) = B(i,j) + Aijk(i,j,k)*mu0(i)^(k-1);
                if k > 1
                    D(i,j) = D(i,j) + Aijk(i,j,k)*(k-1)*mu0(i)^(k-2);
                end
            end
        end
    end
    if isempty(beta0)
        D = [];
        A = [];
        c = [];
    else
        D = diag(D*beta0);
        A = [D -eye(m)];
        c = B*beta0 - nu0;
    end
end
end
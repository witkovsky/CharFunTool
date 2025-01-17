function cf = cfND_FoldedNormal(t, mu, Sigma)
%% cfND_FoldedNormal
%  Computes the characteristic function (CF) of the multivariate folded
%  normal distribution. This distribution arises when the magnitudes (but
%  not the signs) of a multivariate normal random vector are of interest.
%  It is particularly useful in modeling scenarios involving absolute
%  values of normally distributed variables.
%
%  First, let us consider an N-dimensional random vector X that follows a
%  multivariate normal distribution with mean vector mu and positive
%  definite covariance matrix Sigma. Then, the transformation leading to
%  the multivariate folded normal distribution is defined as:
%   W = [|X1|, |X2|, ..., |XN|]',
%  where |Xi| represents the absolute value of the ith component of X.
%
%  The exact form of the characteristic function of the multivariate folded
%  normal distribution is derived in the paper by Benko et al. (2025). For
%  a multivariate normal random vector X with mean vector mu and covariance
%  matrix Sigma, the folded normal random vector W has a characteristic
%  function computed as the sum over all possible sign combinations s.
%  The computation depends on the following components:
%    - s: An N-dimensional vector representing a specific combination of
%      signs (+1 or -1) applied to each dimension of X.
%    - mu_s: The modified mean vector, calculated as mu_s = diag(s) * mu.
%    - Sigma_s: The modified covariance matrix, calculated as Sigma_s =
%      diag(s) * Sigma * diag(s).
%    - Phi_n: The cumulative distribution function (CDF) of the
%      n-dimensional multivariate normal distribution (n = 1, ..., N),
%      which accepts complex arguments and is evaluated over all
%      combinations of the transformed components.
%    - The final computation involves summing the contributions from all
%      sign combinations to accurately account for the folding effect.
%
% SYNTAX:
%  cf = cfND_FoldedNormal(t, mu, Sigma)
%
% INPUTS:
%  t      - (M x N)-matrix of evaluation points in N dimensions.
%  mu     - (N x 1)-vector of means for the normal distribution. Default:
%           [0; ...; 0].
%  Sigma  - (N x N)-covariance matrix. Must be symmetric positive definite.
%           Default: eye(N).
%
% OUTPUT:
%  cf     - (M x 1)-vector of CF values at the specified points in t.
%
% EXAMPLE 1:
% % CF of the bivariate folded normal distribution
%   mu    = [0; 0];
%   Sigma = [1, 0.5; 0.5, 1];
%   szt   = 51;
%   tt    = linspace(-1, 1, szt);
%   [t1, t2] = meshgrid(tt, tt);
%   t     = [t1(:), t2(:)];
%   cf    = cfND_FoldedNormal(t, mu, Sigma);
%   figure; contour(t1, t2, real(reshape(cf, szt, szt)), 'ShowText', 'on');
%   title('Contour Plot of CF of Bivariate Folded Normal Distribution');
%
% EXAMPLE 2:
% % CF of the sum of components of the bivariate folded normal distribution
%   mu    = [0; 0];
%   Sigma = [1, 0.5; 0.5, 1];
%   szt   = 101;
%   t     = linspace(-3, 3, szt);
%   tt    = [t(:), t(:)];
%   cf    = cfND_FoldedNormal(tt, mu, Sigma);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a sum of components of the bivariate folded normal distribution')
%
% EXAMPLE 3:
% % PDF/CDF of the sum of components of the bivariate folded normal distribution
%   mu    = [0; 0];
%   Sigma = [1, 0.5; 0.5, 1];
%   cf    = @(t) cfND_FoldedNormal([t(:), t(:)], mu, Sigma);
%   t_max = 4.4;
%   clear options
%   options.T = t_max;
%   options.xMin = 0;
%   x = linspace(0,30,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 4:
% % CF of a 3-variate folded normal distribution
%   mu    = [1; 2; 3];
%   Sigma = [1, 0.5, 0.1; 0.5, 1, 0.3; 0.1, 0.3, 2];
%   cf    = @(t) cfND_FoldedNormal(t, mu, Sigma);
%   szt   = 5;
%   tt    = linspace(-0.5, 0.5, 11);
%   [t1, t2, t3] = meshgrid(tt);
%   t = [t1(:), t2(:), t3(:)];
%   disp([t, cf(t)]);
%
% EXAMPLE 5:
% % CF of the sum of components of the trivariate folded normal distribution
%   mu    = [1; 2; 3];
%   Sigma = [1, 0.5, 0.1; 0.5, 1, 0.3; 0.1, 0.3, 2];
%   szt   = 51;
%   t     = linspace(-0.3, 0.3, szt);
%   tt    = [t(:), t(:), t(:)];
%   cf    = cfND_FoldedNormal(tt, mu, Sigma);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a sum of components of the trivariate folded normal distribution')
%
% EXAMPLE 6:
% % CF of a 4-variate folded normal distribution
%   mu    = [1; 2; 3; 4];
%   Sigma = [1.0, 0.5, 0.1, 0.0; ... 
%            0.5, 1.0, 0.3, 0.1; ...
%            0.1, 0.3, 1.0, 0.3; ...
%            0.0, 0.1, 0.3, 2.0];
%   cf    = @(t) cfND_FoldedNormal(t, mu, Sigma);
%   szt   = 3;
%   tt    = linspace(-0.2, 0.2, szt);
%   [t1, t2, t3] = meshgrid(tt);
%   t = [t1(:), t2(:), t3(:) t3(:)];
%   disp([t, cf(t)]);
%
% EXAMPLE 7:
% % CF of the sum of components of the 4-variate folded normal distribution
%   mu    = [1; 2; 3; 4];
%   Sigma = [1.0, 0.5, 0.1, 0.0; ... 
%            0.5, 1.0, 0.3, 0.1; ...
%            0.1, 0.3, 1.0, 0.3; ...
%            0.0, 0.1, 0.3, 2.0];
%   szt   = 5;
%   t     = linspace(-0.5, 0.5, szt);
%   tt    = [t(:), t(:), t(:), t(:)];
%   cf    = cfND_FoldedNormal(tt, mu, Sigma);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a sum of components of the 4-variate folded normal distribution')
%
%
% REFERENCES:
% [1] Benko M., Hübnerová Z., Witkovský V. (2025) Characteristic function
%     and moment generating function of multivariate folded normal
%     distribution, Statistical Papers (submitted).
% [2] https://en.wikipedia.org/wiki/Folded_normal_distribution
%
% (c) Zuzana Hubnerova (hubnerova@fme.vutbr.cz)
%     & Matej Benko (benko@fme.vutbr.cz)
% Adapted by
% %   Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: '14-Jan-2025 12:32:25'

%% ALGORITHM
% cf = cfND_FoldedNormal(t,mu,Sigma)

%% Input Validation
narginchk(1, 3);

% Get dimensions of `t`
[M,N] = size(t);

% Set default values for `mu` and `Sigma`
if nargin < 2 || isempty(mu), mu = zeros(N, 1); end
if nargin < 3 || isempty(Sigma), Sigma = eye(N); end

% Validate dimensions
if numel(mu) ~= N
    error('Dimension mismatch: `mu` must be a vector of length N.');
end
if ~isequal(size(Sigma), [N, N])
    error('Dimension mismatch: `Sigma` must be an (N x N) matrix.');
end
if ~issymmetric(Sigma) || any(eig(Sigma) <= 0)
    error('Invalid covariance matrix: `Sigma` must be symmetric positive definite.');
end

%% Determine diagonal block structure of the covariance matrix Sigma

% Create a graph where edges indicate non-zero covariance
g = graph(Sigma ~= 0);

% Find connected components (block structure)
nodeToComponent = conncomp(g);

% Sort nodes by component
[nodeToComponentSorted, p] = sort(nodeToComponent);

%% Compute Characteristic Function

cf = 1;
for id = 1:nodeToComponentSorted(end)
    % indices for the diagonal block of the covariance matrix Sigma
    subp     = p(nodeToComponentSorted == id);
    % Set the selected dimensions from t, mu, and Sigma
    N_id     = length(subp);
    t_id     = t(:,subp);
    mu_id    = mu(subp);
    Sigma_id = Sigma(subp,subp);

    % Pre-compute standard deviations of Sigma_id (unchanged by sign changes)
    % Diagonal elements of Sigma_id define standard deviations
    stdSigma = sqrt(diag(Sigma_id));

    %% Generate all possible combinations of signs (+1, -1)
    % Binary representation of all combinations
    s = dec2bin(0:2^N_id-1) - '0';
    % Replace 0s with -1 to represent sign combinations
    s(s == 0) = -1;

    %% Apply the current sign combination to Sigma and mu
    % Check for equal values to avoid redundant evaluations
    % Initialize a matrix to store information on transformed covariance matrix
    % and mean
    Sigma_mu_check = zeros(2^N_id,N_id^2+N_id);

    for is = 1:2^N_id
        % Set Lambda_s - diagonal matrix with current sign combination
        Lambda_s = diag(s(is, :));
        % Transform Sigma
        Sigma_s = Lambda_s * Sigma_id * Lambda_s;
        % Transform mu
        mu_s = Lambda_s * mu_id;
        % Store the transformed values (flattened Sigma_s and mu_s)
        Sigma_mu_check(is,:) = [transpose(Sigma_s(:)), transpose(mu_s)];
    end

    % Remove duplicate evaluations
    % Use MATLAB's unique function to find unique rows in s_check
    [Sigma_mu_unique, ia, ic] = unique(Sigma_mu_check, 'rows');

    %% Main loop
    % Initialize the vector of characteristic function values cf_id
    cf_id = zeros(M, 1);

    for is = 1:length(ia)
        % Extract mu_s
        mu_s    = Sigma_mu_unique(is, N_id^2+1:end)';        
        % Reshape Sigma_s
        Sigma_s = reshape(Sigma_mu_unique(is, 1:N_id^2), N_id, N_id); 
        % Correlation matrix
        Rho     = Sigma_s ./ (stdSigma * stdSigma');         
        % Compute all d_s = 1i * Sigma_s * tt + mu_s;                   
        d_s     = 1i * transpose(Sigma_s * t_id' + mu_s);
        % Scaled components of d_s      
        scaled_d_s = d_s ./ stdSigma'; 
        % Exponential of the quadratic form plus linear term
        expqf = exp(-0.5 * sum((t_id * Sigma) .* t_id, 2) + 1i * (t_id * mu_s)); 
        % Initial incremental value for the the computed CF
        % Phi is CDF of standard normal distribution that accepts complex arguments 
        value_s = expqf .* Phi(scaled_d_s(:,N_id)); 

        %% Loop over dimensions for multivariate CDF
        for dim = 2:N_id
            % All combinations
            Mnchoosek = nchoosek(1:(N_id-1), dim-1);
            % Diagonal transformation
            D         = -eye(dim);
            D(1, 1)   = 1;
            % Calculate incremental value_s of the CF
            for k = 1:size(Mnchoosek, 1)
                idx = Mnchoosek(k,:);
                value_s = value_s + expqf * (-1)^(dim-1) .* ...
                    nvncdf([scaled_d_s(:,N_id), -scaled_d_s(:,idx)], ...
                    D * Rho([N_id, idx], [N_id, idx]) * D);
            end
        end

        %% Accumulate values of the block characteristic function cf_id
        cf_id = cf_id + sum(ic == ia(is)) * value_s; 
    end

    %%  Combine blocks of the characteristic function
    cf = cf *  cf_id;
end
end

%% Helper function Phi
function p = Phi(z)
% Phi
%  Computes the CDF for the standard normal distribution in complex
%  argument z by using the Faddeeva function.
%
% CDF for the normal distribution.
%  EXAMPLE 1:
%   z = 1i + linspace(-5,5)';
%   f = Phi(z);
%   disp([z f])

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 21-Aug-2023 14:23:58

p = 0.5 * erfcZX(-z / sqrt(2));

end

%% Helper function erfcZX
function f = erfcZX(z)
% erfcZX
%  Computes the complementary error function in complex argument z by using
%  the Faddeeva function.
%

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
%% Faddeeva
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
%% Helper function
function p = nvncdf(b,Rho,tol)
% CDF for the n-dimensional multivariate normal distribution
%  Computes the CDF for the n-dimensional normal distribution in complex
%  argument b and real correlation matrix Rho.
%

% Viktor Witkovsky (witkovsky@gmail.com),
% Zuzana Hubnerova (hubnerova@fme.vutbr.cz)
% Ver.: 26-Mar-2024

if nargin == 2
    tol = 1e-6;
end

%% ALGORITHM
d = size(Rho,1);
if length(b)==d
    b = transpose(b(:));
end
n = size(b,1);

switch d
    case 1
        p = Phi(b);
    case 2
        p = bvncdf(b,Rho(1,2));
    case 3
        p = tvncdf(b,[Rho(1,2),Rho(1,3),Rho(2,3)],tol);
    otherwise

        [~,imin] = min(sum(abs(Rho))); % select a row with the smallest correlations
        if imin > 1 % swap 1 and imin
            id = [imin,2:(imin-1),1,(imin+1):d];
            Rho = Rho(id,id);
            b = b(:,id);
            % imin == 1 correct order
        end

        p1 = Phi(b(:,1)).*nvncdf(b(:,2:end),Rho(2:end,2:end),tol);

        id = find(abs(Rho(1,2:end))>0)+1;
        p2 = zeros(n,length(id),'like',p1);


        for j = id
            loLimit = 0;
            hiLimit = asin(Rho(1,j));
            for in = 1:n % to evaluate for each value in b
                b1 = b(in,1); bj = b(in,j); b2 = b(in,:); b2([1,j]) = [];
                R1j = Rho([1,j],1:end); R1j(:,[1,j]) = [];
                r1j = Rho(1,j);
                A22 = Rho; A22([1,j],:) = []; A22(:,[1,j]) = [];
                if isfinite(b1) && isfinite(bj) && any(~isnan(b2))
                    p2(in,j) = quadgk(@nvnIntegrand,loLimit,hiLimit,'AbsTol',tol/3,'RelTol',0);
                    %else
                    % This piece is zero if either limit is +/- infinity.  If
                    % either is NaN, p1 will already be NaN.
                end
            end
        end


        p = cast(p1 + (sum(p2,2))./(2.*pi), 'like', internal.stats.dominantType(b,Rho));
end
%% nvnIntegrand

    function integrand = nvnIntegrand(theta)
        sintheta = sin(theta);
        cossqtheta = cos(theta).^2; % always positive
        expon = ((b1*sintheta - bj).^2 ./ cossqtheta + b1.^2)/2;

        subnvncdfval = zeros(size(theta));
        for i = 1:length(theta)
            A11inv = [[1 -sintheta(i)];[-sintheta(i) 1]]/cossqtheta(i);
            A1j = [R1j(1,:)*sintheta(i)/r1j;R1j(2,:)] ;
            A = A22 - transpose(A1j)*A11inv*A1j;
            s = sqrt(diag(A));
            subnvncdfval(i) = nvncdf(transpose(transpose(b2) -  transpose(A1j)*A11inv*[b1;bj])./transpose(s),A./(s*transpose(s)));
        end
        integrand = exp(-expon) .* subnvncdfval;
    end

end % nvncdf

%% Helper function
function p = bvncdf(b,rho)
% bvncdf
%  Computes the CDF for the bivariate normal distribution in complex
%  argument b = [z1 z2] and real -1 < rho < 1 (correlation coefficient -
%  here we assume that covariance matrix Sigma is standardized to be a
%  correlation matrix, i.e. Sigma = [1 rho; rho 1]). 
%
%  This is adapted from the Matlab function bvncdf (Copyright 2020-2021 The
%  MathWorks, Inc.) where the original implementation for real arguments is
%  based on Section 2.4 of Genz (2004). 

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 21-Aug-2023 14:23:58

%% ALGORITHM
if length(b)==2 %ZH added to cope with column vectors easily
    b = transpose(b(:));
end
n = size(b,1);
if rho == 0
    p = cast(prod(Phi(b),2), superiorfloat(b,rho));
else
    if abs(rho) < 0.3      % 6 point Gauss Legendre abscissas and weights
        w = [0.4679139345726904  0.3607615730481384  0.1713244923791705];
        y = [0.2386191860831970  0.6612093864662647  0.9324695142031522];
    elseif abs(rho) < 0.75 % 12 point Gauss Legendre abscissas and weights
        w = [0.2491470458134029  0.2334925365383547  0.2031674267230659 ...
             0.1600783285433464  0.1069393259953183  0.04717533638651177];
        y = [0.1252334085114692  0.3678314989981802  0.5873179542866171 ...
             0.7699026741943050  0.9041172563704750  0.9815606342467191];
    else                 % 20 point Gauss Legendre abscissas and weights
        w = [0.1527533871307259  0.1491729864726037  0.1420961093183821  0.1316886384491766  0.1181945319615184 ...
             0.1019301198172404  0.08327674157670475 0.06267204833410906 0.04060142980038694 0.01761400713915212];
        y = [0.07652652113349733 0.2277858511416451  0.3737060887154196  0.5108670019508271  0.6360536807265150 ...
             0.7463319064601508  0.8391169718222188  0.9122344282513259  0.9639719272779138  0.9931285991850949];
    end
    
    if abs(rho) < .925
        p1 = prod(Phi(b),2);
        asinrho = asin(rho);
        w = fliplr([w w]);
        theta = asinrho .* (fliplr([y -y]) + 1) / 2;
        sintheta = sin(theta);
        cossqtheta = cos(theta).^2; % always positive
        h = -b(:,1); k = -b(:,2);
        hk = h .* k;
        ssq = (h.^2 + k.^2) / 2;
        p2 = zeros(size(p1),class(rho));
        for i = 1:n
            if isfinite(hk(i))
                f = exp( -(ssq(i) - hk(i)*sintheta) ./ cossqtheta );
                p2(i) = asinrho * sum(w.*f) / 2;
            else
                % This piece is zero if either limit is +/- infinity.  If
                % either is NaN, p1 will already be NaN.
            end
        end
        p = p1 + p2/(2*pi);

    else % abs(rho) >= .925
        if rho > 0
            p1 = Phi(min(b,[],2));
            p1(any(isnan(b),2)) = NaN;
        else
            p1 = Phi(b(:,1)) - Phi(-b(:,2));
            p1(p1<0) = 0; % max would drop NaNs
        end
        
        s = sign(rho);
        if abs(rho) < 1
            h = -b(:,1); k = -b(:,2);
            shk = s.*h.*k;
            asq = 1 - rho.^2;  a = sqrt(asq);
            b = abs(h - s.*k); bsq = b.^2;
            c = (4 - shk)/8;
            d = c .* (12 - shk)/16;

            t1 = a.*(1 + d.*asq.^2/5 + (c/3 - d.*bsq/15).*(asq - bsq));
            t2 = b.*(1 - c.*bsq/3 + d.*bsq.^2/15);
            p2 = exp(-shk/2 + ...
                log(t1.*exp(-bsq./(2*asq)) - t2.*sqrt(2*pi).*Phi(-b./a)));

            w = [w w];
            x = a .* ([y -y] + 1) / 2; x2 = x.^2; x4 = x2.^2;
            sqrt1mx2 = sqrt(1 - x2);
            t = (sqrt1mx2 - 1) ./ (2*(sqrt1mx2 + 1));
            for i = 1:n
                if isfinite(shk(i))
                    f = exp(-(bsq(i)./x2 + shk(i))/2) ...
                        .* (exp(shk(i).*t)./sqrt1mx2 - (1 + c(i).*x2 + d(i).*x4));
                    p2(i) = p2(i) + a .* sum(w.*f) / 2;
                else
                    % This piece is zero if either limit is +/- infinity.  If
                    % either is NaN, p1 will already be NaN.
                    p2(i) = 0;
                end
            end
        else % abs(rho) == 1
            p2 = zeros(class(rho));
        end
        p = p1 - s.*p2/(2*pi);
    end
end

end % bvncdf
%% Helper function
function p = tvncdf(b,rho,tol)
% CDF for the trivariate normal
% tvncdf
%  Computes the CDF for the trivariate normal distribution in complex
%  argument b = [z1 z2 z3] and real rho = [rho12 rho13 rho23] (correlation
%  coefficients -  here we assume that covariance matrix Sigma is
%  standardized to be a correlation matrix, i.e. Sigma = [1 rho12 rho13;
%  rho21 1 rho23; rho31 rho32 1]). 
%
%  This is adapted from the Matlab function tvncdf (Copyright 2020-2021 The
%  MathWorks, Inc.) where the original implementation for real arguments is
%  based on equation (14) in Section 3.2 of Genz (2004), integrating each
%  term in (14) separately in terms of theta between 0 and asin(rho_j1),
%  using adaptive quadrature.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 21-Aug-2023 14:23:58

if nargin == 2
     tol = 1e-6;
end

%% ALGORITHM
if length(b)==3 %ZH added to cope with column vectors easily
    b = transpose(b(:));
end
n = size(b,1);

% Find a permutation that makes rho_32 == max(rho)
[dum,imax] = max(abs(rho)); %#ok<ASGLU>
if imax == 1 % swap 1 and 3
    rho_21 = rho(3); rho_31 = rho(2); rho_32 = rho(1);
    b = b(:,[3 2 1]);
elseif imax == 2 % swap 1 and 2
    rho_21 = rho(1); rho_31 = rho(3); rho_32 = rho(2);
    b = b(:,[2 1 3]);
else % imax == 3
    rho_21 = rho(1); rho_31 = rho(2); rho_32 = rho(3);
    % b already in correct order
end

p1 = Phi(b(:,1)).*bvncdf(b(:,2:3),rho_32);

if abs(rho_21) > 0
    loLimit = 0;
    hiLimit = asin(rho_21);
    rho_j1 = rho_21;
    rho_k1 = rho_31;
    p2 = zeros(size(p1),'like',p1);
    for i = 1:n
        b1 = b(i,1); bj = b(i,2); bk = b(i,3);
        if isfinite(b1) && isfinite(bj) && ~isnan(bk)
            p2(i) = quadgk(@tvnIntegrand,loLimit,hiLimit,'AbsTol',tol/3,'RelTol',0);
        else
            % This piece is zero if either limit is +/- infinity.  If
            % either is NaN, p1 will already be NaN.
        end
    end
else
    p2 = zeros('like',p1);
end

if abs(rho_31) > 0
    loLimit = 0;
    hiLimit = asin(rho_31);
    rho_j1 = rho_31;
    rho_k1 = rho_21;
    p3 = zeros(size(p1),'like',p1);
    for i = 1:n
        b1 = b(i,1); bj = b(i,3); bk = b(i,2);
        if isfinite(b1) && isfinite(bj) && ~isnan(bk)
            p3(i) = quadgk(@tvnIntegrand,loLimit,hiLimit,'AbsTol',tol/3,'RelTol',0);
        else
            % This piece is zero if either limit is +/- infinity.  If
            % either is NaN, p1 will already be NaN.
        end
    end
else
    p3 = zeros('like',p1);
end

p = cast(p1 + (p2 + p3)./(2.*pi), 'like', internal.stats.dominantType(b,rho));

%% tvnIntegrand

    function integrand = tvnIntegrand(theta)
        % Integrand is exp( -(b1.^2 + bj.^2 - 2*b1*bj*sin(theta))/(2*cos(theta).^2) )
        sintheta = sin(theta);
        cossqtheta = cos(theta).^2; % always positive
        expon = ((b1*sintheta - bj).^2 ./ cossqtheta + b1.^2)/2;

        sinphi = sintheta .* rho_k1 ./ rho_j1;
        numeru = bk.*cossqtheta - b1.*(sinphi - rho_32.*sintheta) ...
            - bj.*(rho_32 - sintheta.*sinphi);
        denomu = sqrt(cossqtheta.*(cossqtheta -sinphi.*sinphi ...
            - rho_32.*(rho_32 - 2.*sintheta.*sinphi)));
        integrand = exp(-expon) .* Phi(numeru./denomu);
    end

end % tvncdf
% RUIN PROBABILITY EXPERIMENTS
% See Kass et al 4.2 The classical ruin process
% Example Example 4.3.7 (Simulating a Poisson ruin process), p. 95
%
% 27-Nov-2016 15:49:45

clear

%% set the CF of the surplus process U(t) = u + c*t - S(t)
% where  where
% U(t) = the insurer’s capital at time t;
% u = U(0) = the initial capital;
% c = the (constant) premium income per unit of time;
% S(t) = X1 + X2 +...+ XN(t),with
% N(t) = the number of claims up to time t, and
% Xi = the size of the ith claim, assumed non-negative.

cfT = @(t) cfX_Exponential(t,1); % CF of the time increments dT
cfX = @(t) cfX_Gamma(t,2,2);     % CF of the individual claims

% CF of the surplus with initial capital u after n claims
% n  number of claims
% u  initial insurer’s capital
% n the (constant) premium income per unit of time
cfU = @(t,n,u,c) exp(1i*t*u) .* cfT(c*t).^n .* cfX(-t).^n;


%% Compute and plot the Ruin probabilities for different number of claims
% Probabilities of negative surplus U (i.e. the the ruin of the insurer)
% after n claims 
N = 300;    % upper limit for number of claims N
RuinProb = zeros(N,1);
RuinLimit = 0;
u  = 12;    % initial insurer’s capital
c  = 1.3;   % the (constant) premium income per unit of time

options.N = 2^10;
options.SixSigmaRule = 6;
for n = 1:N
    cf = @(t)cfU(t,n,u,c);
    [~,ruinProb] = cf2DistGP(cf,RuinLimit,[],options);
    RuinProb(n) = ruinProb;
end

plot(RuinProb,'.')
xlabel('Numer of claims')
ylabel('Probability of Ruin')

%%
% l = 1;
% b = 4;
% psi = @(u) (l/(c*b)) * exp(-(b-l/c)*u);

%%
NonRuinProb = prod(1-RuinProb);
 prob = 0;for i = 1:N, prob = prob + RuinProb(i)*NonRuinProb/(1-RuinProb(i));end
 prob

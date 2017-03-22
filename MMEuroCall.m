function C = MMEuroCall( S, K, r, T, sigma, Q)
%MMEuroCall calculates the current price of a European call option
%
%   MMEuroCall calculates the current price of a European call option
%   based on the assumption that the underlying assett follows a
%   Markov-modulated diffusion process
%
%   C is the current price of the call option
%
%   S is the initial price of the stock
%   K is the strike price
%   r is the risk free rate, in years
%   T is the time until expiry, in years
%   sigma is a vector of states that the volatility can be in (e.g. low and
%    high volatility)
%   Q is the transition rate matrix for the underlying Markov chain
%
% Lachlan Bridges
% (using McKinlay's model from thesis 'Markov-Modulated Models for Derivatives Pricing') 
% 16/01/17


lambda = Q(1,2);
mu = Q(2,1);
v1 = sigma(1);
v2 = sigma(2);
[f1, ~] = OccTime(lambda,mu);

int = @(s) blsprice(S,K,r,T,sqrt(((s*(v1^2)+(T-s)*(v2^2))/T))).*f1(s,T);
C = integral(int,0,T)+exp(-lambda*T)*blsprice(S,K,r,T,v1);




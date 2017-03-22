function [f1, f2] = OccTime(lambda, mu )
%OccTime calculates the conditional (generalised) densities of the times
% spent in the first state of a two state Markov chain, and returns them 
% as functions of x (time spent in state 1) and t (length of the time
% interval)
%
% f1 is the time spent in state 1, given the process in state 1 at time t
% f2 is the time spent in state 1, given the process in state 2 at time t
%
% lambda and mu are the values of the Q transition matrix
%
% Lachlan Bridges
% 09/01/17

f1 = @(x,t) exp(-lambda*x-mu.*(t-x)) .* ( dirac(t-x)+sqrt(lambda*mu.*x...
    ./(t-x)).*besseli(1,2*sqrt(lambda*mu.*x.*(t-x))) +...
    lambda*besseli(0,2*sqrt(lambda*mu.*x.*(t-x))) );

f2 = @(x,t) exp(-lambda*x-mu.*(t-x)) .* ( dirac(x)+sqrt(lambda*mu.*(t-x)...
    ./x).*besseli(1,2*sqrt(lambda*mu.*x.*(t-x))) +...
    lambda*besseli(0,2*sqrt(lambda*mu.*x.*(t-x))) );
end

function f = EmpOccTime(d, Q, y0, T, N, seed)
%EmpOccTime uses the CTMC function to randomly generate Markov chains with
% which it calculates the probability density function of the empirical
% occupation time
%
% f1 is the probability density function, given the process started in
% state y0, function of time
%
% d is the state for which the occupation time shall be found
% Q is the transition rate matrix for the Markov chain
% y0 is the state that the Markov chain started in
% T is the time interval for which the occupation time will be found
% N is the number of Markov chains to generate. larger n -> more accurate
% seed is the seed to be used with the random number generator
%
% Lachlan Bridges
% 16/01/17

[tt,yy] = CTMC(Q,y0,T,N,seed);
tt(tt>T) = T;
tt(tt==0) = T;
tt(isnan(tt)) = T;
tt(:,1) = 0;
tdiffs = tt(:,2:end)-tt(:,1:end-1);

idt = (yy==d);
x = sum(tdiffs.*idt(:,1:(end-1)),2);
f = sort(x);
end


%KSTestOccTime uses the Kolmogorov-Smirnov goodness-of-fit test to test
% whether or not the empirical CDF of the occupation time, found from 
% the Monte Carlo generated Markov chains, match the theoretical
% distribution of the occupation time
%
% Lachlan Bridges
% 02/01/2017

lambda = 10;
nu = 5;
Q=[-lambda, lambda; nu, -nu]; % transition rate matrix

y0 = 1; % intial state
T = 50; % length of time interval
N = 1000; % number of distinct sample paths to generate

x1=EmpOccTime(1,Q,y0,T,N,100);

% function for theoretical pdf of occupation time in state 1
[f1, f2] = OccTime(lambda, nu);
func = @(x) f1(x,T);

x = linspace(0,T,100000); % data points to generate cdf at
x = x(1:end-1); % remove last point since NaN
f = func(x); % theoretical pdf
F = cumsum(f); 
F = F./F(end); % theoretical cdf

cdf = [x' F']; % <- THEORETICAL OCC TIME

[H,P,KSSTAT] = kstest(x1,'CDF',cdf)

% PLOTTING THEORETICAL

s=linspace(0,T,1000);
f1plot = f1(s,T);
f2plot = f1(T-s,T);

plot(s,f1plot, '-b') % state 1 occtime
hold on
plot(s, f2plot, '-r')
title('Probability distribution of occupation times for states 1 and 2')
xlabel('Occupation time')
ylabel('Probability')
%fname = sprintf('../figures/tocctime1_lam%i_nu%i.png',lambda,nu);
%saveas(gcf,fname)
%hold off
figure;

% PLOTTING THEORETICAL VS EMPIRICAL
plot(s,f1plot), hold on
histogram(x1,'Normalization','probability','BinWidth',1)
title(sprintf('Probability distribution of occupation times for state 1\n vs.\n empirical occupation times'))
xlabel('Occupation time')
ylabel('Probability')
%fname = sprintf('../figures/occtime_lam%i_nu%i.png',lambda,nu);
%saveas(gcf,fname)
%hold off
%close all
%MMMilsteinExample tests the MMMilstein function with an example Markov-modulated
% diffusion process with 100 steps per interval
%
% Lachlan Bridges
% 06/01/17

Q=[-5,5;1,-1];
mu = [5, -3];
sigma = [1, 3];

a = @(y,t,J) mu(J);
b = @(y,t,J) sigma(J);
dbdy = @(y,t,J) 0;

x0 = 1;
T = [0,5];
seed = 31;
r = 3; % number of realisations to plot
N=100;
reps=100;

[t,y,J,tt,yy]=MMMilstein(Q,reps,a,b,dbdy,x0,T,N,seed);

plotmmstates(t(1,:),y(1,:),J(1,:),[1,0,0;0,0,1]), hold on
title(sprintf('%Single realisation of a Markov-modulated diffusion\n process, using N=%i time steps per interval',r,reps,N));
xlabel('t','FontSize',14)
ylabel('X_t','FontSize',14)
%saveas(gcf,sprintf('../figures/examplediff.png'));
%close all

figure;
plotmm(t(1:r,:),y(1:r,:)), hold on
[tm, ym] = meanmm(t,y,10000);
plot(tm,ym)

title(sprintf('%i realisations of an MMDP and the mean across\n %i realisations, using N=%i time steps per interval',r,reps,N));
xlabel('t','FontSize',14)
ylabel('X_t','FontSize',14)
%saveas(gcf,sprintf('../figures/examplediffs.png'));
%close all



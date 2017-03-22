function [t, y] = MMEulerMaruyama(Q,reps,a,b,x0,T,N,seed)
%MMEulerMaruyama uses Euler Maruyama method to solve a Markov Modulated SDE
%
% Q is the transition rate matrix of the underlying MC
% reps is number of  simulations to repeat
% a is deterministic term - function of t, y and Z
% b is random term - function of t, y and Z
% x0 is starting pos
% T is final time
% N is number of discrete time intervals to use
% seed is rng seed
%
% t is the 1xN vector of time steps that are used
% y is the (reps)xN vector of solutions at time steps
%
% Lachlan Bridges
% 06/01/17

randn('seed',seed)

dt=T/N;
t=0:dt:T;

[tt, yy]=CTMC(Q,1,T,reps,seed);
for j=1:reps
    for i=1:length(t)
        J(j,i)=yy(j,find(tt(j,:)>t(i),1)-1);
    end
end

y=zeros(reps,N+1);

y(:,1)=repmat(x0,reps,1);
Z=randn(reps,N);
for i=2:N+1
    y(:,i)=y(:,i-1)+a(y(:,i-1),i*dt,J(i))*dt + b(y(:,i-1),i*dt,J(i))*sqrt(dt).*Z(:,i-1);
end

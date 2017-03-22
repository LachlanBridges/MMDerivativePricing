function [t, y] = EulerMaruyama(reps,a,b,x0,T,N,seed)
%EulerMaruyama uses the Euler Maruyama method to solve an SDE
%
% reps is number of  simulations to repeat
% a is deterministic term - function of t and y
% b is random term - function of t and y
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

y=zeros(reps,N+1);
dt=T/N;
t=0:dt:T;

y(:,1)=repmat(x0,reps,1);
Z=randn(reps,N);
for i=2:N+1
    y(:,i)=y(:,i-1)+a(y(:,i-1),i*dt)*dt + b(y(:,i-1),i*dt)*sqrt(dt).*Z(:,i-1);
end

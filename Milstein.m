function [t, y] = Milstein(reps,a,b,dbdy,x0,T,N,seed)
%Milstein uses the Milstein method to solve an SDE
%
% reps is number of  simulations to repeat
% a is deterministic term - function of t and y
% b is random term - function of t and y
% dbdy is the derivative of b w.r.t y - function of t and y
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

rng(seed)

T0 = T(1);
Tf = T(2);
y=zeros(reps,N+1);
dt=(Tf-T0)/N;
t=T0:dt:Tf;

y(:,1)=repmat(x0,reps,1);
Z=randn(reps,N);
for i=2:N+1
    y(:,i)=y(:,i-1)+a(y(:,i-1),i*dt+T0)*dt + b(y(:,i-1),i*dt+T0)*sqrt(dt).*Z(:,i-1)...
        +dbdy(y(:,i-1),i*dt+T0).*b(y(:,i-1),i*dt+T0).*(Z(:,i-1).^2-1)*dt/2;
end

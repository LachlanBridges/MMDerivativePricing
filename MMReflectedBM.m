function [t, X, B] = MMReflectedBM(Q,reps,a,b,x0,T,N,seed)
%MMReflectedBM generates reflected brownian motions along with its
%non-reflected counterpart.
%
% t is the 1xN vector of time steps that are used
% X is the (reps)xN vector of solutions  to the REFLECTED SDE at t
% B is the (reps)xN vector of solutions  to the NON-REFLECTED SDE at t
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
% Credit to:  Hanbook of Monte Carlo Methods (Kroese, Taimre, Botev)
% Adapted by: Lachlan Bridges
% 06/01/17
 
randn('seed',seed)
dt = T/N;
t=0:dt:T;

X=zeros(reps,N+1); M=X; B=X;

% initial conditions
B(:,1)=x0;
X(:,1)=x0;

% generate CTMC paths
[tt, yy]=CTMC(Q,1,T,reps,seed);
for j=1:reps
    for i=1:length(t)
        J(j,i)=yy(j,find(tt(j,:)>t(i),1)-1);
    end
end

for r=1:reps
    for k=2:N+1
        Y=b(X(r,k-1),t(k-1),J(r,k))*sqrt(dt)*randn;
        U=rand(1);
        B(r,k) = B(r,k-1) + a(X(r,k-1),t(k-1),J(r,k))*dt -Y;
        M = (Y + sqrt(Y^2 - 2*dt*log(U)))/2;
        X(r,k) = max(M-Y, X(r,k-1)+dt*a(X(r,k-1),t(k-1),J(r,k))-Y);
    end
end

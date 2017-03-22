function [tt, yy] = CTMC(Q,y0,T,reps, seed)
%CTMC generates sample paths for a Continuous Time Markov Chain
%
% tt is the vector of transition times
% yy is the vector of states
%
% Q is the transition rate matrix of the Markov Chain
% y0 is the initial state of the Markov Chain
% T is the length of time to generate
% reps is the number of sample paths to generate
% seed is the rng seed
%
% Credit to:  Hanbook of Monte Carlo Methods (Kroese, Taimre, Botev)
% 06/01/17


rng(seed)
S=size(Q,2);
q=-diag(Q);
K=diag(1./q)*Q + eye(S);

for j=1:reps
    n=1;
    t=0; y=y0;
    yy(j,1)=y;
    tt(j,1)=t;
    while t < T
        A = -log(rand)/q(y);
        y = find(cumsum(K(y,:))>rand, 1 );
        t = t + A;
        n=n+1;
        tt(j,n) = t;
        yy(j,n) = y;
    end
end
yy(tt==0) = nan;
tt(tt==0) = nan;
yy(:,1) = y0;
tt(:,1) = 0;
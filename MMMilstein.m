function [t, y, J, tt, yy] = MMMilstein(Q,reps,a,b,dbdy,x0,T,N,seed)
%MMMilstein uses the Milstein method to solve a Markov Modulated SDE
%
% Q is the transition rate matrix of the underlying MC
% reps is number of  simulations to repeat
% a is deterministic term - function of t, y and Z
% b is random term - function of t, y and Z
% dbdy is the derivative of b w.r.t y - function of t,y,Z
% x0 is starting pos
% T is a 2 element vector, specifying the start and end time
% N is number of discrete time intervals to use between each transition
% seed is rng seed
%
% t is the (reps)xN vector of time steps that are used
% y is the (reps)xN vector of solutions at time steps
%
% Lachlan Bridges
% 06/01/17

rng(seed)

[tt, yy]=CTMC(Q,1,T(2),reps,seed);
maxjumps = (size(yy,2));
intervals = (maxjumps-1)*(N);
t = zeros(reps,intervals);
y = t;
J = y;
for i = 1:reps
    Jrep = [];
    trep = [];
    yrep = [];
    ttt = tt(i,:);
    yyy = yy(i,:);
    jumps = find(isnan(yyy),1);
    if isempty(jumps);
        jumps = maxjumps;
    end
    for j = 1:jumps-1
        newT = [ttt(j), ttt(j+1)];
        state = yyy(j);
        if j==1
            y0 = x0;
        else
            y0 = yrep(end);
        end
        [tempt,tempy]=Milstein(1,@(y,t) a(y,t,state),@(y,t) b(y,t,state),...
            @(y,t) dbdy(y,t,state),y0,newT,N,seed);
        tempt(1) = [];
        tempy(1) = [];
        tempJ = repmat(state,size(tempt));
        trep = [trep tempt];
        yrep = [yrep tempy];
        Jrep = [Jrep tempJ];
    end
    trep = [trep zeros(1,intervals-length(trep))];
    yrep = [yrep zeros(1,intervals-length(yrep))];
    Jrep = [Jrep zeros(1,intervals-length(Jrep))];
    t(i,:) = trep;
    y(i,:) = yrep;
    J(i,:) = Jrep;
end
y(t==0 | t>T(2)) = nan;
J(t==0 | t>T(2)) = nan;
t(t==0 | t>T(2)) = nan;
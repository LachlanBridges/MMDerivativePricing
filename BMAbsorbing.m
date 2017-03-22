function X = BMAbsorbing(B,lb,ub)
%BMAbsorbing generates a diffusion process with absorbing boundary 
% conditions, given a normal diffusion process as an input
%
% X is the (reps)xN vector of solutions to the diffusion process with jump
%
% B is the original diffusion process
% lb is the value of the lower boundary at which the process will be
% absorbed (put [] if there is no lower bound)
% ub is the value of the upper boundary at which the process will be
% absorbed (put [] if there is no lower bound)
% 
% Lachlan Bridges
% 06/01/17

% case when no upper/lower boundary
if isempty(lb)
    lb = -inf;
end
if isempty(ub)
    ub = inf;
end

N = size(B,2);
reps = size(B,1);

X=B;

for r=1:reps
    for k=2:N
        if X(r,k) >= ub
            X(r,k:end) = ub;
        end
        if X(r,k) <= lb
            X(r,k:end) = lb;
        end
    end
end
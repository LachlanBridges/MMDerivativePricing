function X = BMInstantJump(B,lb,ljump,ub,ujump)
%BMInstantJump generates a diffusion process with instantaneous jump
%boundary conditions, given a normal diffusion process as input
%
% X is the (reps)xN vector of solutions to the diffusion process with jump
%
% B is the original diffusion process
% lb is the value of the lower boundary at which the jump will occur
% ljump is the value that is jumped to upon hitting the lower boundary
% ub is the value of the upper boundary at which the jump will occur
% ujump is the value that is jumped to upon hitting the upper boundary
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

ujs = ub - ujump;
ljs = lb - ljump;

for r=1:reps
    for k=2:N
        if X(r,k) >= ub
            X(r,k:end) = X(r,k:end)-ujs;
        end
        if X(r,k) <= lb
            X(r,k:end) = X(r,k:end)-ljs;
        end
    end
end
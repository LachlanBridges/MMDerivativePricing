function X = BMReflect(B,lb)
%BMReflect generates reflected brownian motions given a non-reflected
%brownian motion as input
%
% X is the (reps)xN vector of solutions  to the REFLECTED SDE at t
%
% B is the NON-REFLECTED SDE
% lb is the value of the lower bound for the REFLECTED SDE
%
% Lachlan Bridges
% 06/01/17
 
N = size(B,2);
reps = size(B,1);

X=zeros(reps,N);

%initial condition
X(:,1)=B(1,1);


for r=1:reps
    for k=2:N
        infm = -1*(min(B(r,1:k))-lb);
        X(r,k) = B(r,k) + max(0,infm);
    end
end
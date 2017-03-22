function [] = plotmmstates(t,y,J,cols)
%plotmmstates plots each realisation of a markov modulated diffusion with
% different color for the diffusion in each state

% t is the input of times. each row represents a separate realisation
% y is the input of values. each row represents a separate realisation
% J is the input of states. each row represents a separate realisation
% cols is a (states)x3 vector of colors for the plot. columns represent rgb

if nargin<4
    rng('shuffle')
    num = size(t,1);
    states = max(J(:));
    cols = rand(states,3);
else
    num = size(t,1);
    states = max(J(:));
end

hold on
for i = 1:num
    tt = t(i,:);
    yy = y(i,:);
    JJ = J(i,:);
    for j = 1:states
        ttj = tt;
        yyj = yy;
        ttj(~(JJ==j))=nan;
        yyj(~(JJ==j))=nan;
        plot(ttj,yyj,'color',cols(j,:))
    end
end
hold off
end


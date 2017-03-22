function [] = plotmm(t,y,line)
%plotmm plots each realisation of a markov modulated diffusion in a
% different color on a plot

% t is the input of times. each row represents a separate realisation
% y is the input of values. each row represents a separate realisation
% line is optional input specifying color and type of line

if nargin < 3
  line = [];
end

num = size(t,1);
hold on
for i = 1:num
    if isempty(line)
        plot(t(i,:),y(i,:))
    else
        plot(t(i,:),y(i,:),line)
    end
end
hold off
end


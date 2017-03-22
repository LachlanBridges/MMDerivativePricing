function [ tf, yf] = meanmm(t, y, N)
%meanmm finds the mean across multiple realisations of a diffusion process
% using interpolation

% t is a matrix of times, rows represents realisations
% y is a matrix of values, rows represent realisations
% N is the number of total timesteps

yf = [];
tf = linspace(min(t(:)),max(t(:)),N);
for i=1:size(t,1)
    tt = t(i,:);
    yy = y(i,:);
    yy(isnan(tt)) = [];
    tt(isnan(tt)) = [];
    yi = interp1(tt,yy,tf);
    yf = [yf ; yi];
end
yf = mean(yf);
end


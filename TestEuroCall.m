%TestEuroCall uses our implementation (MMEuroCall) of McKinlay's model and
% compares the resulting price with the prices calculated by McKinlay to
% ensure our model is correct
%
% Lachlan Bridges
% 16/01/17

S0 = 50; % intial stock price
r = 0.1; % risk free interest rate
sigma = [0.25, 0.55]; % volatility
T = 1; % expiry

%strike prices to test:
K = [29.943175, 39.952443, 49.828298, 59.899056, 69.829658, 79.921838];

% correct prices (from McKinlay's thesis):
m1 = [23.229596, 15.354187, 9.3150276, 5.2848860, 2.9956081, 1.7274812];
m2 = [23.568903, 16.353291, 10.919362, 7.0734063, 4.5946202, 2.9876057];
m3 = [23.104774, 15.020033, 8.8040807, 4.7166456, 2.4739770, 1.3000675];
oprices = [m1 ; m2 ; m3];

% model paramaters:
lambda = [1, 3, 1];
nu = [1, 1, 3];

for i=1:3
    fprintf('MODEL %i:\n\n',i)
    Q = [-lambda(i), lambda(i); nu(i), -nu(i)]; % transition matrix
    for j=1:length(K)
        C=MMEuroCall(S0,K(j),r,T,sigma,Q);
        fprintf('McKinlay''s price: %f,     Our Price: %f,     Difference:%f\n',...
            oprices(i,j),C,abs(C-oprices(i,j)))
    end
    fprintf('\n')
end



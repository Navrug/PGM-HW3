data = importdata('EMGaussian.data');
testdata = importdata('EMGaussian.test');

k = 4;
rng(0);

pi0 = 0.25 * ones(1, 4);
A = ones(4,4)/6 + diag(ones(4,1))/3;

t = 10;
T = 500;

eps = 1e-3;
[ ~, ~, mu, sigma, ~ ] = GM2(data,4,eps);

%%
p = compute_p(alphas(data, A, mu, sigma, pi0), betas(data, A, mu, sigma));
[~,labels] = max(p,[],2);

% Plot the state probabilities.
figure
for i = 1:k
   subplot(k, 1, i)
   plot(1:100, p(1:100,i))
end

% Plot the geometrical positions with most probable state.
% This is not the same thing as the most probable global state.
figure
scatter(data(:,1), data(:,2), [], labels, 'x'); hold on;
scatter(mu(:,1), mu(:,2), [], 'black', 'filled');  hold on;

%%
viter = viterbi( data, A, mu, sigma, pi0 );
%mean(viterbi == labels') % = 71 percent

% Plot the state probabilities.
figure
for i = 1:k
   subplot(k, 1, i)
   plot(1:100, p(1:100,i))
end

% Plot the geometrical positions with most probable state.
% This is not the same thing as the most probable global state.
figure
scatter(data(:,1), data(:,2), [], viter', 'x'); hold on;
scatter(mu(:,1), mu(:,2), [], 'black', 'filled');  hold on;

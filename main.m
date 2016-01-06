train = importdata('EMGaussian.data');
test = importdata('EMGaussian.test');

k = 4;
T = 500;

%% Question 2
pi = 0.25 * ones(1, 4);
A = ones(4,4)/6 + diag(ones(4,1))/3;

% Use parameters from HW2.
eps = 1e-3;
rng(0);
[ ~, ~, mu, sigma] = GM2(train, NaN, 4, eps);

% Forward-backward algorithm.
a = alphas(train, A, mu, sigma, pi);
b = betas(train, A, mu, sigma);
pz = compute_pz(a, b);

% Plot the state probabilities.
figure
for i = 1:k
   subplot(k, 1, i)
   plot(1:100, pz(1:100,i))
end

%% Question 4 : EM
rng(0);
eps = 1e-3;
[labels, pi, A, mu, sigma, l_train, l_test, GM_train, GM_test] = EM(train, test, k, eps);

% Plot the geometrical positions.
figure
scatter(train(:,1), train(:,2), [], labels, 'x'); hold on;
scatter(mu(:,1), mu(:,2), [], 'black', 'filled');  hold on;

% Plot the train and test likelihoods.
figure
plot(1:length(l_train), l_train, 1:length(l_test), l_test); hold on;

% Compare final likelihoods of GM and HMM.
fprintf('GM has likelihoods %f and %f\n', GM_train, GM_test)
fprintf('HMM has likelihoods %f and %f\n', l_train(length(l_train)), l_test(length(l_test)))

%% Question 8 : Viterbi

viter = viterbi(train, A, mu, sigma, pi);

% Plot the state probabilities.
% --> it is not an output of the viterbi algo (only states)

% Plot the geometrical positions with most probable state.
% This is not the same thing as the most probable global state.
figure
scatter(train(:,1), train(:,2), [], viter', 'x'); hold on;
scatter(mu(:,1), mu(:,2), [], 'black', 'filled');  hold on;

%% Question 9 

% % Forward-backward algorithm.
a = alphas(test, A, mu, sigma, pi);
b = betas(test, A, mu, sigma);
pz = compute_pz(a, b);

% Plot the state probabilities.
figure
for i = 1:k
   subplot(k, 1, i)
   plot(1:100, pz(1:100,i))
end
hold off;

%% Question 10

% Plot the states
[~, states] = max(pz, [], 2);
figure
scatter(1:100, states(1:100))

%% Question 11

viter_test = viterbi(test, A, mu, sigma, pi);

% Plot the viterbi states along with most likely states according to the
% marginal probability
hold on;
scatter(1:100, viter_test(1:100), '+'); hold off;


%% Question 12

max_k = 5;
likelihoods = zeros(max_k);
for k = 2:max_k
    rng(0);
    eps = 1e-3;
    [labels, pi, A, mu, sigma, l_train, ~, ~, ~] = EM(train, NaN, k, eps);
    likelihoods(k) = l_train(length(l_train));
end
figure
plot(1:max_k, likelihoods)
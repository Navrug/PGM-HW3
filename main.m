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
a = alpha( t, data, k, A, mu, sigma, pi0 );
b = beta(  t, T, data, k, A, mu, sigma );
function [ q ] = viterbi( data, A, mu, sigma, pi )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

k = length(A);
T = length(data);

logv = zeros(T,k);
x = zeros(T, k);

logv(1,:) = log(pi) + log(mvnpdf(data(1,:),mu,sigma)');

for t=2:length(data)
   for i=1:k
       [logv(t,i), x(t,i)] = max( log(A(:,i)') + logv(t-1,:) );
       logv(t,i) = logv(t,i) + log(mvnpdf(data(t,:),mu(i,:),sigma(:,:,i)));
   end
end

q = zeros(1, T);
[~, q(T)] = max(x(T,:));
for t=T:-1:2
    q(t-1)=x(t,q(t));
end

end


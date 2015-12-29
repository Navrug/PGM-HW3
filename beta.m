function [ beta ] = beta(  t, T, data, k, A, mu, sigma )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

p_cond = zeros(1, k);
for i=1:k
    p_cond(i) = mvnpdf(data(t,:),mu(i,:),sigma(:,:,i));
end

betas = zeros(T, k);
betas(T, :) = ones(1, k);

for s=1:(T-t)
    for i=1:k
        p_cond(i) = mvnpdf(data(T-s,:),mu(i,:),sigma(:,:,i));
    end
    betas(T-s,:) = (betas(T-s+1,:)*A).*p_cond;
end
beta = betas(t,:);



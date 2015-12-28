function [ alpha ] = alpha( t, data, k, A, mu, sigma, pi )
% pi is a line

p_cond = zeros(1, k);
for i=1:k
    p_cond(i) = mvnpdf(data(t,:),mu(i,:),sigma(:,:,i));
end

alphas = zeros(t,4);
alphas(1,:) = pi .* p_cond;

for s=2:t
    alphas(s,:) = (alphas(s-1,:)*A).*p_cond;
end

alpha = alphas(t,:);


function [ likelihood ] = cl( data, labels, mu, sigma )

[n, ~] = size(data);
likelihood = 0;
for i=1:n
    likelihood = likelihood - log(mvnpdf(data(i,:), mu(labels(i),:), sigma(:,:,labels(i))));
end

end


function [ labels, pi, mu, sigma ] = GM2( data, k , eps )
%GM2 Summary of this function goes here
%   Detailed explanation goes here

   [labels, mu] = k_means(data,k);
   
   [n,d] = size(data);
   pi = zeros(1,k);
   sigma = zeros(d,d,k);
   for i = 1:k
      pi(i) = sum(labels==i)/n;
      sigma(:,:,i) = eye(d);
   end
   
   L = -1e10;
   T = 100;
   t = 0;
   while t < T
      t = t + 1;  
      
      % E step
      L_old = L;
      L = 0;
      tau = zeros(n,k);
      for i = 1:n
         for j = 1:k
            tau(i,j) = pi(j)*mvnpdf(data(i,:), mu(j,:), sigma(:,:,j));
         end
         L = L + log(tau(i,labels(i)));
         tau(i,:) = tau(i,:)/sum(tau(i,:));
      end
     
      if L - L_old < eps
         break
      end
      
      % M step
      pi = mean(tau,1);
      mu = zeros(k,d);
      for j = 1:k
         mu(j,:) =  + tau(:,j)'*data/sum(tau(:,j))  ; 
         sigma(:,:,j) = zeros(d);
         for i = 1:n 
            sigma(:,:,j) = sigma(:,:,j) + tau(i,j)*(data(i,:)-mu(j,:))'*(data(i,:)-mu(j,:));
         end
         sigma(:,:,j) = sigma(:,:,j)/(sum(tau(:,j)));
      end
      [~,labels] = max(tau,[],2);
   end
   
end


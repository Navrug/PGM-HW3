function [ labels, pi, mu, sigma, L, L_test ] = GM2( data, test, k , eps )
%GM2 Summary of this function goes here
%   Detailed explanation goes here
   L_test = NaN;
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
         L = L + log(sum(tau(i,:)));
         tau(i,:) = tau(i,:)/sum(tau(i,:));
      end
      
      
      if ~isnan(test)
          L_test = 0;
          test_tau = zeros(n,k);
          for i = 1:n
             for j = 1:k
                test_tau(i,j) = pi(j)*mvnpdf(test(i,:), mu(j,:), sigma(:,:,j));
             end
             L_test = L_test + log(sum(test_tau(i,labels(i))));
          end
      end
      fprintf('L = %f\n', L)
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
   fprintf('GM EM performed %i iterations\n', t)
   
end


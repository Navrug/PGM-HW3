function [ labels, pi, A, mu, sigma ] = EM( data, k, eps )
%   Summary of this function goes here
%   Detailed explanation goes here

   [~, pi, mu, sigma] = GM2(data, k, eps);
   A = ones(k,k)/6 + diag(ones(k,1))/3;

   [T,d] = size(data);

   L = -1e10;
   max_it = 100;
   it = 0;
   while it < max_it
      it = it + 1;
      % E step
      loga = alphas(data,A,mu,sigma,pi);
      logb = betas(data,A,mu,sigma);
      pz = compute_pz(loga, logb);
      pzz = compute_pzz(loga, logb, data, A, mu, sigma);
      
      L_old = L;
      L = pz(1,:)*pi';
      for t = 2:T
         L = L + pz(t,:)*log(mvnpdf(data(t,:),mu,sigma));
         L = L + sum(sum(A.*reshape(pzz(t-1,:,:),k,k)));
      end
      fprintf('L = %f\n', L)
      if L - L_old < eps
         break
      end
      
      % M step
      A = reshape(sum(pzz, 1),k,k);
      A = diag(1 ./ sum(A,2))*A;
      pi = mean(pz,1);
      mu = zeros(k,d);
      for j = 1:k
         mu(j,:) =  + pz(:,j)'*data/sum(pz(:,j))  ; 
         sigma(:,:,j) = zeros(d);
         for t = 1:T 
            sigma(:,:,j) = sigma(:,:,j) + pz(t,j)*(data(t,:)-mu(j,:))'*(data(t,:)-mu(j,:));
         end
         sigma(:,:,j) = sigma(:,:,j)/(sum(pz(:,j)));
      end
   end
   fprintf('EM performed %i iterations.\n', it)
   [~,labels] = max(pz,[],2);
   
end


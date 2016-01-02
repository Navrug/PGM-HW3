function [ labels, pi, A, mu, sigma, train_likelihoods, test_likelihoods, GM_train, GM_test ] = EM( train, test, k, eps )
%   Summary of this function goes here
%   Detailed explanation goes here

   [~, pi, mu, sigma, GM_train, GM_test] = GM2(train, test, k, eps);
   A = ones(k,k)/6 + diag(ones(k,1))/3;

   [T,d] = size(train);
   
   train_likelihoods = [];
   if ~isnan(test)
    test_likelihoods = [];
   end

   L = -1e10;
   max_it = 100;
   it = 0;
   while it < max_it
      it = it + 1;
      % E step
      loga = alphas(train,A,mu,sigma,pi);
      logb = betas(train,A,mu,sigma);
      pz = compute_pz(loga, logb);
      pzz = compute_pzz(loga, logb, train, A, mu, sigma);
      
      L_old = L;
      L = pz(1,:)*pi';
      for t = 2:T
         L = L + pz(t,:)*log(mvnpdf(train(t,:),mu,sigma));
         L = L + sum(sum(A.*reshape(pzz(t-1,:,:),k,k)));
      end
      train_likelihoods = horzcat(train_likelihoods, L);
      fprintf('L = %f\n', L)
      
      if ~isnan(test)
          test_loga = alphas(test,A,mu,sigma,pi);
          test_logb = betas(test,A,mu,sigma);
          test_pz = compute_pz(test_loga, test_logb);
          test_pzz = compute_pzz(test_loga, test_logb, test, A, mu, sigma);

          test_L = test_pz(1,:)*pi';
          for t = 2:T
             test_L = test_L + test_pz(t,:)*log(mvnpdf(test(t,:),mu,sigma));
             test_L = test_L + sum(sum(A.*reshape(test_pzz(t-1,:,:),k,k)));
          end
          test_likelihoods = horzcat(test_likelihoods, test_L);
          fprintf('L_test = %f\n', test_L)
      end
      
      if L - L_old < eps
         break
      end
      
      % M step
      A = reshape(sum(pzz, 1),k,k);
      A = diag(1 ./ sum(A,2))*A;
      pi = mean(pz,1);
      mu = zeros(k,d);
      for j = 1:k
         mu(j,:) =  + pz(:,j)'*train/sum(pz(:,j))  ; 
         sigma(:,:,j) = zeros(d);
         for t = 1:T 
            sigma(:,:,j) = sigma(:,:,j) + pz(t,j)*(train(t,:)-mu(j,:))'*(train(t,:)-mu(j,:));
         end
         sigma(:,:,j) = sigma(:,:,j)/(sum(pz(:,j)));
      end
   end
   fprintf('EM performed %i iterations.\n', it)
   [~,labels] = max(pz,[],2);
   
end


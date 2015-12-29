function [ p ] = compute_pzz( loga, logb, data, A, mu, sigma )
%P Summary of this function goes here
%   Detailed explanation goes here
   [~,k] = size(loga);
   
   logp = zeros(length(loga)-1, k, k);
   % p(:,i,j) = A(i,j)
   for t=1:(length(loga)-1)
      logp(t,:,:) = log(A);
   end
   
   for j = 1:k
      % p(t,i,j) *= b(t+1,j)
      logp(:,:,j) = logp(:,:,j) + logb(2:length(logb),:);
      for t=1:(length(loga)-1)
         % p(t,i,j) *= p(i|j)
         logp(t,:,j) = log(mvnpdf(data(t,:),mu,sigma)');
      end
   end
   sz = size(loga);
   for i = 1:k
      % p(t,i,j) *= a(t,i)
      logp(:,i,:) = logp(:,i,:) + reshape(loga(1:(length(loga)-1),:), sz(1)-1, 1, sz(2));
   end
   % sum(t) = p(t,1,1)
   logsum = logp(:,1,1);
   for i = 1:k
      for j = (1+(i==1)):k
         % sum(t) += p(t,i,j)
         logsum = logaddexp(logsum, logp(:,i,j));
      end
   end
   % p(t,i,j) /= sum(t)
   p = exp(bsxfun(@minus,logp, logsum));
end


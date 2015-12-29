function [ p ] = compute_pzz( loga, logb, data, A, mu, sigma )
%P Summary of this function goes here
%   Detailed explanation goes here
   [~,k] = size(loga);
   
   logp = zeros(length(loga)-1, k, k);
   for t=1:(length(loga)-1)
      logp(t,:,:) = log(A);
   end
   for j = 1:k
      logp(:,:,j) = logp(:,:,j) + logb(2:length(logb),:);
      for t=1:(length(loga)-1)
         logp(t,:,j) = log(mvnpdf(data(t,:),mu,sigma)');
      end
   end
   sz = size(loga);
   for i = 1:k
      logp(:,i,:) = logp(:,i,:) + reshape(loga(1:(length(loga)-1),:), sz(1)-1, 1, sz(2));
   end
   logsum = logp(:,1,1);
   for i = 1:k
      for j = (1+(i==1)):k
         logsum = logaddexp(logsum, logp(:,i,j));
      end
   end
   p = exp(bsxfun(@minus,logp, logsum));
end


function [ p ] = compute_pzz( loga, logb, data, A, mu, sigma )
%P Summary of this function goes here
%   p(zt, zt+1|y) = a(t)p(yt+1|zt+1)b(t+1)A(qt, qt+1) / normalisation

   [~,k] = size(loga);
   
   % Computing unnormalised probability.
   logp = zeros(length(loga)-1, k, k);
   % p(:,zt,zt+1) = A(zt,ztp1)
   for t=1:(length(loga)-1)
      logp(t,:,:) = log(A);
   end
   for zt = 1:k
      for t=1:(length(loga)-1)
         % p(t,zt,ztp1) *= p(ytp1|ztp1)
         logp(t,zt,:) = logp(t,zt,:) + reshape(log(mvnpdf(data(t+1,:),mu,sigma)), 1,1,k);
      end
   end
   sz = size(logb);
   % p(t,zt,ztp1) *= b(t+1,ztp1)
   logp = bsxfun(@plus, logp, reshape(logb(2:length(logb),:), sz(1)-1, 1, sz(2)));
   % p(t,i,j) *= a(t,i)
   logp = bsxfun(@plus, logp, reshape(loga(1:(length(loga)-1),:), sz(1)-1, sz(2),1));
   
   % Normalisation.
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


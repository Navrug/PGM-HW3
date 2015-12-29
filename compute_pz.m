function [ p ] = compute_pz( loga, logb )
%P Summary of this function goes here
%   Detailed explanation goes here
   [~,k] = size(loga);
   % sum = a1*b1
   logsum = loga(:,1) + logb(:,1);
   for j=2:k
      % sum += aj*bj
      logsum = logaddexp(logsum, loga(:,j)+logb(:,j));
   end
   % p = a*b/sum(a*b)
   p = exp(bsxfun(@minus,loga + logb, logsum));
end


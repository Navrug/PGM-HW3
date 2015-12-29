function [ p ] = compute_pz( loga, logb )
%P Summary of this function goes here
%   Detailed explanation goes here
   [~,k] = size(loga);
   logsum = loga(:,1) + logb(:,1);
   for j=2:k
      logsum = logaddexp(logsum, loga(:,j)+logb(:,j));
   end
   p = exp(bsxfun(@minus,loga + logb, logsum));
end


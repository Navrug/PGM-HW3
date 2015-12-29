function [ p ] = compute_p( loga, logb )
%P Summary of this function goes here
%   Detailed explanation goes here
   [~,k] = size(loga);
   logsum = loga(:,1) + logb(:,1);
   for j=2:k
      logsum = logsum + log(1 + exp(loga(:,j) + logb(:,j) - logsum));
   end
   p = exp(bsxfun(@minus,loga + logb, logsum));
end


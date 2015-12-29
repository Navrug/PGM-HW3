% function [ betas ] = betas( data, A, mu, sigma, pi )
% %BETAS Summary of this function goes here
% %   Detailed explanation goes here
% 
% betas = zeros(length(data), length(A));
% betas(length(data),:) = 1;
% for t=(length(data)-1):-1:1
%    p_yz = mvnpdf(data(t+1,:),mu,sigma)'; % = p(y|z)
%    betas(t,:) = (p_yz.*betas(t+1,:))*A';
% end


function [ logb ] = betas( data, A, mu, sigma )
% Outputs all the alpha messages, a txk message.
k = length(A);
logb = zeros(length(data), k);
for t=(length(data)-1):-1:1
   for i=1:k
      logb(t,i) = log(A(i,1))+logb(t+1,1);
      for j=2:k
         logb(t,i) = logb(t,i) + log(1 + exp(A(i,j)+logb(t+1,j) - logb(t,i)));
      end
   end
   logp_yz = log(mvnpdf(data(t+1,:),mu,sigma)'); % = p(y|z)
   logb(t,:) = logp_yz + logb(t,:);
end
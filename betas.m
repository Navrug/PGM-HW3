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
   logp_yz = log(mvnpdf(data(t+1,:),mu,sigma)'); % = p(y|z)
   for i=1:k
      % b(t,i) = A(i,1)*b(t+1,1)
      logb(t,i) = log(A(i,1))+logb(t+1,1)+logp_yz(1);
      for j=2:k
         % b(t,i) = A(i,j)*b(t+1,j)*p_yz(j)
         logb(t,i) = logaddexp(logb(t,i), log(A(i,j))+logb(t+1,j)+logp_yz(j));
      end
   end
end
% function [ alphas ] = alphas( data, A, mu, sigma, pi )
% % Outputs all the alpha messages, a txk message.
% 
% alphas = zeros(length(data), length(A));
% alphas(1,:) = pi.*mvnpdf(data(1,:),mu,sigma)';
% for t=2:length(data)
%    p_yz = mvnpdf(data(t,:),mu,sigma)'; % = p(y|z)
%    alphas(t,:) = p_yz.*(alphas(t-1,:)*A);
% end

function [ loga ] = alphas( data, A, mu, sigma, pi )
% Outputs all the alpha messages, a txk message.
k = length(A);
loga = zeros(length(data), k);
% a(1,i) = pi*p(y1|i)
loga(1,:) = log(pi) + log(mvnpdf(data(1,:),mu,sigma)');
for t=2:length(data)
   for i=1:k
      % a(t,i) = A(1,i)*a(t-1,1)
      loga(t,i) = log(A(1,i))+loga(t-1,1);
      for j=2:k
         % a(t,i) += A(j,i)*a(t-1,j)
         loga(t,i) = logaddexp(loga(t,i), log(A(j,i))+loga(t-1,j));
      end
   end 
   logp_yz = log(mvnpdf(data(t,:),mu,sigma)'); % = p(y|z)
   % a = p_yz * a
   loga(t,:) = logp_yz + loga(t,:);
end
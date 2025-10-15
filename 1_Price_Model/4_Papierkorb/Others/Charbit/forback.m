function [alfa,beta,ct,LLr] = forback(N,O,Tr,A,p,B)

%FORBACK Forward-Backward Procedure
%   [alfa beta ct scale] = forback(N,O,A,Sigma,Mu,p,B)
%   N     - number of states
%   O     - observation sequence (T x d x R)
%   Tr    - length of each occurence (1 x R)
%   A     - state transition matrix (N x N)
%   p     - initial state distrib. (N x 1)
%   B     - log-likelyhood of observations (T x N x R)
%   --- output ---
%   alfa  - forward variable (T x N x R)
%   beta  - backward variable (T x N x R)
%   ct    - scaling factor of each occurence (T x R)
%   LLr   - log likelihood (R x 1) 
% NO other function is used
%
%Last changed 01 09 99

[T d R] = size(O);

alfa = zeros(T,N,R);
beta = zeros(T,N,R);
ct = zeros(T,R);
LLr = zeros(R,1);

for r=1:R
% initialization utile pour un calcul de la log-vrais
   alfa(1,:,r) = p .* B(1,:,r); 
   ct(1,r) = sum(alfa(1,:,r));      
   alfa(1,:,r) = alfa(1,:,r) / ct(1,r);
   for t=2:Tr(r)
      for s=1:N
         alfa(t,s,r)=(alfa(t-1,:,r)*A(:,s))*B(t,s,r);
      end
      ct(t,r) = sum(alfa(t,:,r));             % determine scaling factor
      alfa(t,:,r) = alfa(t,:,r) / ct(t,r);    % apply scaling
   end
end
%== pour les calculs ultérieurs de gamma et xi, on peut initialiser 
% \beta_{T,1:N,r} avec n'importe quelle valeur
for r=1:R,
   beta(Tr(r),:,r) = ones(1,N);
end
for r=1:R,
   for t=Tr(r)-1:-1:1,
      for s=1:N
         beta(t,s,r)=(beta(t+1,:,r).*B(t+1,:,r)) * A(s,:)';
      end
      beta(t,:,r)=beta(t,:,r) /ct(t+1,r);
   end
end
for r=1:R,
   LLr(r) = -sum(log(ct(1:Tr(r),r))) ; 
end
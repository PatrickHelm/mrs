function [gamma,sum_t_gamma] = gammafc(alfa,beta,Tr)

%GAMMAFC Update the a priori probability
%   [gamma sum_t_gamma] = gammafc(alfa,beta,Tr)
%   alfa        - forward variable (T x N x R)
%   beta        - backward variable (T x N x R)
%   Tr          - length of each occurence (1 x R)
%   --- output ---
%   gamma       - state probability given the output sequence. 
%                 Size = T x N x R
%   sum_t_gamma - sum on time of alfa*beta. Size N x R
%
% NO other function is used
%
%Last changed 09 09 99

[T N R] = size(alfa);
gamma=zeros(T,N,R);
sum_t_gamma = zeros(N,R);

for r=1:R,
   for tr=1:Tr(r),
      gamma(tr,:,r) = alfa(tr,:,r) .* beta(tr,:,r);
      sum_t_gamma(:,r) = sum_t_gamma(:,r) + gamma(tr,:,r)';
   end
end
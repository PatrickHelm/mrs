function mu_next = mufc(O,Tr,gamma,gsum)

%MUFC Updates the mean matrix (T x d)
%   mu_next = mufc(O,gamma,gsum)
%   O     - observation sequence (T x d x R)
%   Tr    - length of each occurence (1 x R)
%   gamma - a priori probability
%   gsum  - sum of gamma
% NO other function is used
%
%Last changed 01 09 99

[T d R] = size(O);
[N R]=size(gsum);

mu_next = zeros(N,d);

for s=1:N
   for u=1:d
      gs=0;
      for r=1:R,
         mu_next(s,u) = mu_next(s,u) + gamma(1:Tr(r),s,r)'*O(1:Tr(r),u,r);
         gs = gs + gsum(s,r);
      end
      mu_next(s,u)=mu_next(s,u)/gs;
   end
end
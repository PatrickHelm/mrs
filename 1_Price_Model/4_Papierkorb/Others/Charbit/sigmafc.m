function sigma_next = sigmafc(O,Tr,gamma,gsum,Mu)

%SIGMAFC Updates covariance matrix
%   sigma_next = sigmafc(O,gamma,gsum,mu)
%   O     - observation sequence (T x d x R)
%   Tr    - length of the non null observation
%   gamma - a priori probability (T x N x R)
%   gsum  - sum of gamma (N x R)
%   Mu    - mean matrix (N x d)
%   --- output ---
%   sigma_next - variance (N x d)
%
% NO other function is used
%
%Last changed 01 09 99

[T d R] = size(O);
[N d] = size(Mu);

sigma_next = zeros(N,d);

for n=1:N,
   for k=1:d,
      sigma_next(n,k)=0;
      for r=1:R,
         for t=1:Tr(r),
            Xtirn = (O(t,k,r)-Mu(n,k));
            sigma_next(n,k) = sigma_next(n,k)+gamma(t,n,r)*Xtirn*Xtirn;
         end
      end
   end
   sigma_next(n,:)=sigma_next(n,:) / sum(gsum(n,:));
end
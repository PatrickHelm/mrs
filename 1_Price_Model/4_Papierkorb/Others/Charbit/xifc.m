function xi_next = xifc(A,B,alfa,beta,Tr,ct)

%XIFC Updates xi
%   xi_next = xifc(A,B,alfa,beta,N,Tr)
% ---- input ------
%   A       - state transition probability distrib (N x N)
%   B       - observation probability given state (T x N x R)
%   alfa    - forward variable (T x N x R)
%   beta    - backward variable (T x N x R)
%   Tr      - length of each occurence (1 x T)
% ---- output -----
%   xi_next - joint state probability given observation. 
%             Size = (T-1) x N x N x R 
% NO other function is used
%
%Last changed 09 09 99

[T N R] = size(alfa);
xi_next = zeros(T-1,N,N,R);

for r=1:R,
   for tr=1:Tr(r)-1,
      for sf=1:N
         for si=1:N
            xi_next(tr,si,sf,r)=alfa(tr,si,r)*beta(tr+1,sf,r)*B(tr+1,sf,r)*A(si,sf);
         end
      end
      xi_next(tr,:,:,r)=xi_next(tr,:,:,r)/ct(tr+1,r);
   end
end
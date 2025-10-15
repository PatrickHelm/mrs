function aij = aijfc(xi,N,R)

%AIJFC Updates the state transition matrix
%   aij = aijfc(xi,N)
%
%   xi     - joint probability given observation (T-1) x N x N x R
%   N      - number of states
%   R      - number of occurences
%
%  aij     - transition matrix (N x N)
%
% NO other function is used
%
%Last changed 02 09 99

aij=zeros(N,N);
for si=1:N
   for sf=1:N
      for r=1:R
         aij(si,sf)=aij(si,sf)+sum(xi(:,si,sf,r));
      end
   end
   aij(si,:)=aij(si,:)/sum(aij(si,:));
end
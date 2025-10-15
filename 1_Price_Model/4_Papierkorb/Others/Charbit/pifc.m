function pi_next = pifc(gsum)

%PIFC Updates initial state probability
%   gsum - sum of a priori probability
%
%   pi_next = pifc(gsum,T,R)
% NO others fucntions are used
%
%Last changed 01 09 99

[N R]=size(gsum);
pi_next=zeros(1,N);

for k=1:N,
   pi_next(k) = sum(gsum(k,:));
end
pi_next=pi_next/sum(pi_next);
function [Pi,A,Mu,Sigma] = hmminit(O,Tr,N)

%HMMINIT Initializes the variables for the HMM reestimation procedure.
%   [Pi A Mu Sigma] = hmminit(O,Tr,N)
%   O     - ONE observation sequence (T x d)
%   Tr    - length of the non null observation (length of each observation)
%   N     - number of states

%   --- output ---
%   Pi    - inital state distribution (1 x N)
%   A     - state transition matrix (N x N)
%   Mu    - mean vector (N x d) (mu of the regimes)
%   Sigma - variance (N x d) (sigma of the regimes)

%   --- input ---
O=[1,2,3,3,2,3,3,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4] %Price observations
Tr=1 %length of the non null observation 
N=2 %Number of states

%   --- start ---
[T d] = size(O); %length and dimension of observation sequence
Otmp = O(1:Tr,:);
% cut Tr in N blocks
TronN=fix(Tr/N);
for s=1:N %for all states
   Ox=O((s-1)*TronN+1:s*TronN,:);
   Mu(s,:)=mean(Ox); %Mu for all regimes
   Sigma(s,:)=std(Ox) .^2; %Mu for all regimes
end

Pi = zeros(1,N);
Pi(1) = 1;		%always start in state 1
A = zeros(N,N); %state transition matrix
for i=1:N,
   A(i,i:N) = ones(1,N-i+1)/(N-i+1);
end                
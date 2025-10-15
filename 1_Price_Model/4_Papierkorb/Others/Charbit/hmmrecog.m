function LL = hmmrecog(O, Tr, Mu, Sigsq, A, Pi)

%HMMRECOG Recognition of an observation sequence
%   LL = hmmrecog(O,Tr,Mu,Sigsq,A,Pi)
%   O      - observation sequence to compare to the model
%   Tr     - length of each observation
%   --- model parameters ---
%   Mu     - mean matrix
%   Sigsq  - variance
%   A      - state transition matrix
%   Pi     - initial state distribution
%   --- output ---
%   LL     - log likelihood of the observation
% Used functions
%          - GAUSSIAN
%          - FORBACK

%Input: 
O=[1,2,3,3,2,3,3,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4] %Price observations
Tr=1 %length of the non null observation 
N=2 %Number of states

% --- Declaration ---
[T d R] = size(O);
[N d] = size(Mu);

alfa = zeros(T,N,R);
beta = zeros(T,N,R);
gamma = zeros(T,N,R);

% compute gaussian distribution
B = gaussian(Mu,Sigsq,O,Tr);

% update forward-backward
[alfa beta ct LL] = forback(N,O,Tr,A,Pi,B);
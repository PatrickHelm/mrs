%

function [Mu,Sigma,A,Pi,LL,ct] = hmmtrain(O,Tr,N,noRep,tol)

% HMMTRAIN Training of the HMM for multiple data sets.
%   [Mu,Sigma,A,Pi,LL,ct] = hmmtrain(O,Tr,N,noRep,tol)
%   O      - R observation sequences of cepstral coef. (T x d x R)
%   where       d      - number of cepstral coeffs (typical : 12)
%               T      - max on R of Tr(r=1:R)
%               R      - number of occurences of the word to be learned
%   Tr     - length of each occurence (1 x R)
%   N      - number of states
%   --- optional parameters ---
%   noRep  - number of repetitions for the reestimation
%   tol    - tolerance (default 0.001)
%   --- output ---
%   Mu      - mean matrix (N x d)
%   Sigma   - covariance matrix (N x d)
%   A       - state transition probability distrib (N x N)
%   B       - observation likelyhood function (T x d x R)
%   LL      - loglikelihood of the R joint occurences at
%             each step
%   ct      - scaling factor

% Used functions
%           - HMMINIT
%           - GAUSSIAN
%           - FORBACK
%           - GAMMAFC
%           - XIFC
%           - MUFC
%           - SIGMAFC
%           - PIFC
%           - AIJFC

%   --- input ---
O=[1,2,3,3,2,3,3,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4] %Price observations
Tr=1 %length of the non null observation 
N=2 %Number of states
noRep=1
tol=0.001

if nargin<5   tol=0.01; end;
if nargin<4   noRep=30; end;

% --- Declaration ---
[T d R] = size(O);
alfa = zeros(T,N,R);
beta = zeros(T,N,R);
gamma = zeros(T,N,R);
gammasum=zeros(N,R);
LL = zeros(noRep,1);
xi = zeros(T-1,N,N,R);
B  = zeros(T,N,R);
Mu  = zeros(N,d);
Sigma = zeros(N,d);
lgv = -inf;
upd_pi = 0; % if 0 we don't re-estime Pinitial=[1 0 .. 0]

% --- Initialization ---
disp('Initialization');
[Pi A Mu Sigma] = hmminit(O(:,:,1),Tr(1),N);
%
% --- Reestimation ---
disp('Starting reestimation');

for nr=1:noRep,

   % compute gaussian distribution
   B = gaussian(Mu,Sigma,O,Tr);
   
   % update forward-backward
   [alfa beta ct LLr] = forback(N,O,Tr,A,Pi,B);

   % update gamma
   [gamma gammasum] = gammafc(alfa,beta,Tr);
   
   % state probability
   xi = xifc(A,B,alfa,beta,Tr,ct);
    
   % update mean
   Mu = mufc(O,Tr,gamma,gammasum);
  
   % update covariance
   Sigma = sigmafc(O,Tr,gamma,gammasum,Mu);
      
   % update initial state probabilities
   if upd_pi,
      Pi = pifc(gammasum); %
   end
   
   % update state transition matrix   
   A = aijfc(xi,N,R);
  
   % determine Log Likelihood
   lgv_p = lgv;
   lgv = sum(LLr);
   disp(sprintf('Iteration no: %2i, Likelyhood = %2.6f',nr, lgv));
   LL(nr)=lgv;
   if(abs(lgv-lgv_p)<tol)	% if tolerance reached: exit iteration
      break;
   end
end
LL=LL(1:nr);
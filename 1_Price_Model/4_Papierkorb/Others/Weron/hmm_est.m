% Janczura: HMM_EST: Matlab function to estimate parameters of a 2-state
% Hidden markov Model (HMM)
function [Ksi_tT,Param,P,Ksi_t1t_10,LogL]=hmm_est(Data,Model,startParam,startP,startKsi_t1t_10)
        % startP: initial transition matrix
        % P: trained transition matrix
        % Param: estimated parameters
        % Ksi_tT: returned smoothed inferences
        % Ksi_t1t_10: probabilities classifying the first observation to
        % one of the regimes 
        % LOGL: log-likelihood of the fitted model
        % Model: specifies distribution in each regime:
                % Model='G-G': Gaussian distribution in both regimes
                % Model='LN-LN': Lognormal in both regimes
                % Model='G-LN': Gaussian in first, lognormal in second r.
        % Data: Time series data
        
%   The first column (KSI_TT) or row (PARAM, P) contains results for the
%   first regime and the second column/row for the second regime.
%   [KSI_TT,PARAM,P,KSI_T1T_10,LOGL]=HMM_EST(DATA,MODEL) additionally 
%   returns probabilities KSI_T1T_10 classifying the first observation to
%   one of the regimes and log-likelihood LOGL of the fitted model.
%
%   HMM_EST(DATA,MODEL,STARTPARAM,STARTP,STARTKSIT1T_10) 
%   allows to specify initial (model-dependent) parameter estimates 
%   STARTPARAM (a 2x5 vector), initial transition matrix STARTP and initial 
%   estimates STARTKSIT1T_10 of probabilities classifying the first 
%   observation. Default values for the latter three parameters are:
%       switch Model
%         case 'G-G  '
%             startParam = [NaN, 1, .1, NaN, NaN; NaN, 3, .5, NaN, NaN];
%         case 'LN-LN'
%             startParam = [NaN, 1, .5, 0, NaN; NaN, 3, 1, 0, NaN];
%         case 'G-LN '
%             startParam = [NaN, 1, .5, 0, NaN; NaN, 3, 1, 0, NaN];   
%      end
%       STARTP = [0.8, 0.2; 0.6, 0.4];
%       STARTKSIT1T_10 = [0.6, 0.4];
%
%   See also MRS2IR_EST, PS2R_EST
%

% Input:
   Data=[1,2,3,4,3,2,3,4,5,2,3,4,3];
   Model='G-G'
   startParam = [NaN, 1, .1, NaN, NaN; NaN, 3, .5, NaN, NaN];
   startP = [0.8, 0.2; 0.6, 0.4];
   startKsi_t1t_10 = [0.6, 0.4];
   
if nargin < 5
    startKsi_t1t_10 = [.6, .4]; % ???
end
if nargin < 4
    startP = [.8, .2; .6, .4]; %initial transition probability matrix
end
if nargin <3
    startParam=nan*ones(2,5);
	switch Model
        case 'G-G  '
            startParam = [nan, 1, .1, nan, nan; nan, 3, .5, nan, nan];
        case 'LN-LN'
            startParam = [nan,1,.5, 0, nan; nan, 3, 1, 0, nan];
        case 'G-LN '
            startParam = [nan,1,.5, 0, nan; nan, 3, 1, 0, nan]; 
    end     
end

% Initial estimation step
[newKsi_tT, newKsi_tt, newP, newKsi_t1t_10, newLogL] = mrs_EM_Smoother(Data, Model, startParam, startP, startKsi_t1t_10);   
newParam = mrs_EM_MLE(Data, Model, newKsi_tT);

% Estimation loop
iteration = 1;
dif = 1;
while and(dif > 10^(-8), iteration < 100)
    oldP = newP;
    oldKsi_t1t_10 = newKsi_t1t_10;
    oldParam = newParam;
    oldKsi_tT = newKsi_tT;
    oldLogL = newLogL;
    [newKsi_tT, newKsi_tt, newP, newKsi_t1t_10, newLogL] = mrs_EM_Smoother(Data, Model, oldParam, oldP, oldKsi_t1t_10);    
    newParam = mrs_EM_MLE(Data, Model, newKsi_tT);
    
    diffP = max(abs(newP-oldP));
    diffKsi_t1t_10 = max(abs(newKsi_t1t_10-oldKsi_t1t_10));
    diffParam = max(abs(newParam-oldParam));
    diffKsi_tT = max(abs(newKsi_tT-oldKsi_tT));
    diffLogL = abs(newLogL-oldLogL);
    dif = max([diffP diffKsi_t1t_10 diffParam diffKsi_tT diffLogL]);
    iteration = iteration+1;
end

Ksi_tT = newKsi_tT;
Param = newParam;
P = newP;
ind=find(P<0);
if ~isempty(ind)
    P(ind)=0;
end
Ksi_t1t_10 = newKsi_t1t_10;
LogL = newLogL;

%=========================================================================
% mrs_EM_MLE
%=========================================================================

function newParam = mrs_EM_MLE(Data, Model, Ksi_tT)
% MRS_EM_MLE ML estimation of regime processes parameters.

sData = Data(2:end);
sKsi_tT = Ksi_tT(2:end, 2);
bData = Data(1:end);
bKsi_tT = Ksi_tT(1:end, 1);
switch Model
    case 'G-G  '
        sNewC = sum((sData).*sKsi_tT) / sum(sKsi_tT);
        sNewSigma2 = sum(((sData)-sNewC).^2 .* sKsi_tT) / sum(sKsi_tT);
        sNewParam = [nan, sNewC, sNewSigma2, nan, nan];
        bNewC = sum((bData).*bKsi_tT) / sum(bKsi_tT);
        bNewSigma2 = sum(((bData)-bNewC).^2 .* bKsi_tT) / sum(bKsi_tT);
        bNewParam = [nan, bNewC, bNewSigma2, nan, nan];
    case 'G-LN '
        qS=0;
        ind=find(sData>qS);
        sNewC = sum(log(sData(ind)-qS).*sKsi_tT(ind)) / sum(sKsi_tT(ind));
        sNewSigma2 = sum((log(sData(ind)-qS)-sNewC).^2 .* sKsi_tT(ind)) /sum(sKsi_tT(ind));
        sNewParam = [nan,sNewC, sNewSigma2, qS, nan];
        bNewC = sum((bData).*bKsi_tT) / sum(bKsi_tT);
        bNewSigma2 = sum(((bData)-bNewC).^2 .* bKsi_tT) / sum(bKsi_tT);
        bNewParam = [nan, bNewC, bNewSigma2, nan, nan];
    case 'LN-LN'
        qS=0;
        ind=find(sData>qS);
        sNewC = sum(log(sData(ind)-qS).*sKsi_tT(ind)) / sum(sKsi_tT(ind));
        sNewSigma2 = sum((log(sData(ind)-qS)-sNewC).^2 .* sKsi_tT(ind)) /sum(sKsi_tT(ind));
        sNewParam = [nan,sNewC, sNewSigma2, qS, nan];
        ind=find(bData>qS);
        bNewC = sum(log(bData(ind)-qS).*bKsi_tT(ind)) / sum(bKsi_tT(ind));
        bNewSigma2 = sum((log(bData(ind)-qS)-bNewC).^2 .* bKsi_tT(ind)) /sum(bKsi_tT(ind));
        bNewParam = [nan,bNewC, bNewSigma2, qS, nan];
end
newParam = [bNewParam; sNewParam];

%==========================================================================
% mrs_EM_Smoother 
%==========================================================================

function [Ksi_tT, Ksi_tt, newP, newKsi_t1t_10, LogL] = mrs_EM_Smoother(Data, Model, oldParam, oldP, oldKsi_t1t_10)  

N = length(Data);
Ksi_tt = [oldKsi_t1t_10; zeros(N-1, 2)];
Ksi_t1t = [oldKsi_t1t_10; zeros(N-1, 2)];
Ksi_tT = [oldKsi_t1t_10; zeros(N-1, 2)];

switch Model
    case 'LN-LN'
        qS=0;
        sOldC = oldParam(2,2);
        sOldSigma2 = oldParam(2,3);
        bOldC = oldParam(1,2);
        bOldSigma2 = oldParam(1,3);
        Eta = [nan,nan;lognSpdf(Data(2:end), bOldC, sqrt(bOldSigma2),qS) , lognSpdf(Data(2:end), sOldC, sqrt(sOldSigma2),qS)]; 
        for row = 2:N
            Ksi_tt(row, :) = Eta(row, :) .* Ksi_t1t(row-1, :) / sum( Eta(row, :) .* Ksi_t1t(row-1, :) );
            Ksi_t1t(row, :) = Ksi_tt(row, :) * oldP;
        end
    case 'G-G  '
        sOldC = oldParam(2,2);
        sOldSigma2 = oldParam(2,3);
        bOldC = oldParam(1,2);
        bOldSigma2 = oldParam(1,3);
        Eta = [nan,nan;normpdf(Data(2:end), bOldC, sqrt(bOldSigma2)), normpdf(Data(2:end), sOldC, sqrt(sOldSigma2))]; 
        for row = 2:N
            Ksi_tt(row, :) = Eta(row, :) .* Ksi_t1t(row-1, :) / sum( Eta(row, :) .* Ksi_t1t(row-1, :) );
            Ksi_t1t(row, :) = Ksi_tt(row, :) * oldP;
        end
    case 'G-LN '
        qS=0;
        sOldC = oldParam(2,2);
        sOldSigma2 = oldParam(2,3);
        bOldC = oldParam(1,2);
        bOldSigma2 = oldParam(1,3);
        Eta = [nan,nan;normpdf(Data(2:end), bOldC, sqrt(bOldSigma2)), lognSpdf(Data(2:end), sOldC, sqrt(sOldSigma2),qS)]; 
        for row = 2:N
            Ksi_tt(row, :) = Eta(row, :) .* Ksi_t1t(row-1, :) / sum( Eta(row, :) .* Ksi_t1t(row-1, :) );
            Ksi_t1t(row, :) = Ksi_tt(row, :) * oldP;
        end
end
Ksi_tT(end, :) = Ksi_tt(end, :);

for row = (N-1):-1:2
    Ksi_tT(row, :) = Ksi_tt(row, :) .* ( ( Ksi_tT(row+1, :) ./ Ksi_t1t(row, :) ) * oldP' );
end
%Log-likelihood
Aux = Eta(2:end,:) .* Ksi_tT(2:end,:);
LogL = sum(log(Aux(:,1) + Aux(:,2)));

%Transition probabilities
p_11 = sum( oldP(1,1) * Ksi_tT(3:end, 1) .* Ksi_tt(2:end-1, 1) ./ Ksi_t1t(2:end-1, 1) ) / sum( Ksi_tT(2:end-1, 1) );
p_22 = sum( oldP(2,2) * Ksi_tT(3:end, 2) .* Ksi_tt(2:end-1, 2) ./ Ksi_t1t(2:end-1, 2) ) / sum( Ksi_tT(2:end-1, 2) );
newP = [p_11, 1-p_11; 1-p_22, p_22];

newKsi_t1t_10 =  Ksi_tT(2,:);

%==========================================================================
% mrs_EM_Smoother auxiliary functions 
%==========================================================================

function Y = lognSpdf(X, mu, sigma,c)
Y = zeros(length(X),1);
ind = find(X>c);
Y(ind) = lognpdf_mod(X(ind)-c,mu,sigma);

function Y = lognpdf_mod(X, mu, sigma)
Y = lognpdf(X, mu, sigma);
for i=1:length(X)
    if isnan(Y(i))
        Y(i) = 0;
    end
end

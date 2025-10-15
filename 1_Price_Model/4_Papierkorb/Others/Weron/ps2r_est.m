% Janczura, Weron: Matlab function to estimate parameters of a 2-regime
% parameter switching (PS) model

function [Ksi_tT,Param,P,Ksi_t1t_10,LogL]=ps2r_est(Data,gammas,display,startParam,startP,startKsi_t1t_10)
%PS2R_EST Estimates parameters of a 2-regime parameter switching (PS) model. 
%   [KSI_TT,PARAM,P]=PS2R_EST(DATA) returns smoothed inferences KSI_TT, 
%   estimated parameters PARAM and transition matrix P of a 2-regime 
%   parameter switching (PS) model, i.e. a 2-regime Markov regime-switching 
%   (MRS) model with both regimes driven by AR(1) processes of the form: 
%   X(t+1)=phi_i*X(t)+c_i+sigma_i*|X(t)|^g_i*N(0,1).
%   The first column (KSI_TT) or row (PARAM, P) contains results for the
%   base regime and the second column/row for the spike regime.
%   [KSI_TT,PARAM,P,KSI_T1T_10,LOGL]=PS2R_EST(DATA) additionally returns 
%   probabilities KSI_T1T_10 classifying the first observation to one of 
%   the regimes and log-likelihood LOGL of the fitted model. 
%   
%   [...]=PS2R_EST(DATA,GAMMAS) allows to specify the known values of g_1
%   and g_2 (as 1x2 vector [g_1, g_2]; default value GAMMAS=[]). If vector 
%   GAMMAS is not given (or GAMMAS=[]) g_1 and g_2 are estimated, otherwise 
%   a numerically less demanding algorithm is used for estimating the 
%   remaining parameters.
%
%   PS2R_EST(DATA,GAMMAS,DISPLAY,STARTPARAM,STARTP,STARTKSIT1T_10) allows 
%   to specify initial (model-dependent) parameter estimates STARTPARAM 
%   (a 2x5 vector), initial transition matrix STARTP and initial estimates
%   STARTKSIT1T_10 of probabilities classifying the first observation.
%   Default values for the latter three parameters are:
%       STARTPARAM = [0.3, 15, 1, 0, NaN; 0.3, 15, 1, 0, NaN];
%       STARTP = [0.8, 0.2; 0.6, 0.4];
%       STARTKSIT1T_10 = [0.6, 0.4];
%   The third pameter, DISPLAY (default DISPLAY=1), is a flag which defines 
%   whether the calibration results are dispayed in the command window (1) 
%   or not (0).
%   
%   Example:
%       Param = [0.2,2,1,0;0.4,3,1,1]; P=[0.5,0.5;0.4,0.6];
%       [Y,S] = ps2r_sim(P,Param,3.5,[1,0],1000,1);
%       ps2r_est(Y);
%       ps2r_est(Y,[0,1]);  
%   
%   See also PS2R_SIM
%
%   Reference(s):
%   [1] J.Janczura, R.Weron (2011) Efficient estimation of Markov 
%   regime-switching models: An application to electricity spot prices. 
%   Working paper version available at: 
%   http://ideas.repec.org/p/wuu/wpaper/hsc1102.html

%   Written by Joanna Janczura and Rafal Weron (2011.02.28)
%   Revised by Joanna Janczura and Rafal Weron (2011.10.03)

if nargin < 6
    startKsi_t1t_10 = [0.6, 0.4];
end
if nargin < 5 
    startP = [0.5, 0.5; 0.2, 0.8];
end
if nargin<2
    gammas = [];
end
if nargin < 4
    startParam = [0.3, 15, 1, 0, NaN; 0.3, 15, 1, 0, NaN];
    T = mean(Data);
    s = std(diff(Data));
    Pind = (Data>(T+s));
    startParam = mrs_EM_MLE(Data, [1-Pind,Pind], startParam, gammas);
end
if nargin<3
    display = 1;
end
if numel(gammas)==2,
    % Assign known values of g_1 and g_2
    startParam(1,4) = gammas(1);
    startParam(2,4) = gammas(2);
else
    % Estimate g_1 and g_2
    gammas = [];
end

% Initial estimation step 
[newKsi_tT,  newP, newKsi_t1t_10, newLogL] = mrs_EM_Smoother(Data, startParam, startP, startKsi_t1t_10);   
newParam = mrs_EM_MLE(Data, newKsi_tT, startParam, gammas);

% Estimation loop 
iteration = 1;
dif = 1;
while and(dif > 10^(-8), iteration < 100)
    oldP = newP;
    oldKsi_t1t_10 = newKsi_t1t_10;
    oldParam = newParam;
    oldKsi_tT = newKsi_tT;
    oldLogL = newLogL; 
    [newKsi_tT, newP, newKsi_t1t_10, newLogL] = mrs_EM_Smoother(Data, oldParam, oldP, oldKsi_t1t_10);    
    newParam = mrs_EM_MLE(Data, newKsi_tT, oldParam, gammas);
    
    diffP = max(abs(newP-oldP));
    diffKsi_t1t_10 = max(abs(newKsi_t1t_10-oldKsi_t1t_10));
    diffParam = max(abs(newParam-oldParam));
    diffKsi_tT = max(abs(newKsi_tT-oldKsi_tT));
    diffLogL = abs(newLogL-oldLogL);
    dif = max([diffP diffKsi_t1t_10 diffParam diffKsi_tT diffLogL]);
    iteration = iteration+1;
end

% Results
Ksi_tT = newKsi_tT;
Param = newParam;
P = newP;
ind = find(P<0);
if ~isempty(ind)
    P(ind) = 0;
end
Ksi_t1t_10 = newKsi_t1t_10;
LogL = newLogL;
% Display results in the command window
if display
    Summary = mrs_Summary(Param, P, LogL);
    disp([' ']);
    disp(['Gamma=' sprintf('%0.2g',Param(1,4)) ' for base regime']);
    disp(['Gamma=' sprintf('%0.2g',Param(2,4)) ' for spike regime']);
    % mean-reverting process
    disp([' ']);        
    disp(['Two state regime switching model with mean-reverting process for spikes, LogL=' num2str(LogL)]);                        
    disp([' ']);        
    disp(['regime   phi_i   c_i      sigma^2_i  E(Y_{t,i})  Var(Y_{t,i})  q_ii     P(R=i)']);
    disp(['base    ' sprintf('%7.5f',Summary(1,1)) '  ' sprintf('%7.5f',Summary(1,2)) '  ' sprintf('%7.5f',Summary(1,3)) '    ' ...
        sprintf('%7.5f',Summary(1,4)) '     ' sprintf('%7.5f',Summary(1,5)) '       ' sprintf('%7.5f',Summary(1,6)) '  ' sprintf('%7.5f',Summary(1,7))]);
    disp(['spike   ' sprintf('%7.5f',Summary(2,1)) '  ' sprintf('%7.5f',Summary(2,2)) '  ' sprintf('%7.5f',Summary(2,3)) '    ' ...
        sprintf('%7.5f',Summary(2,4)) '     ' sprintf('%7.5f',Summary(2,5)) '       ' sprintf('%7.5f',Summary(2,6)) '  ' sprintf('%7.5f',Summary(2,7))]);    
end

%=========================================================================
% Internally used routine(s)
%=========================================================================

%=========================================================================
% mrs_Summary
%=========================================================================

function Summary = mrs_Summary(Param,P,logL)

% Mean and variance
bParam = Param(1,:);
bPhi = bParam(1);
bC = bParam(2);
bSigma2 = bParam(3);
bMean =  bC / (1-bPhi);
bVar = nan;
sParam = Param(2,:);
sPhi = sParam(1);
sC = sParam(2);
sSigma2 = sParam(3);
sMean = sC / (1-sPhi);
sVar = nan;
% Unconditional probabilities
p_11 = P(1,1);
p_22 = P(2,2);
P_Rt_1 = (1-p_22) / (2-p_11-p_22);
P_Rt_2 = (1-p_11) / (2-p_11-p_22);
% Summary matrix
% Mean-reverting process
sPhi = sParam(1);
sC = sParam(2);
sSigma2 = sParam(3);        
Summary = [ bPhi,   bC,     bSigma2,  bMean,  bVar,   p_11,   P_Rt_1;
            sPhi,   sC,     sSigma2,  sMean,  sVar,   p_22,   P_Rt_2;
            logL,     nan,    nan,      nan,    nan,    nan,    nan];
                
%=========================================================================
% mrs_EM_MLE
%=========================================================================

function newParam = mrs_EM_MLE(Data,Ksi_tT,oldParam,gammas)
% MRS_EM_MLE ML estimation of regime processes parameters. 
sData = Data(1:end);
sKsi_tT = Ksi_tT(2:end, 2);

bData = Data(1:end);
bKsi_tT = Ksi_tT(2:end, 1);

if isempty(gammas),
    % Estimate g_1 and g_2
    % Spike regime
    sNewG = fminsearch(@(g) mle_g(sData,sKsi_tT,g),oldParam(2,4),optimset('Display','off'));
    sNewParam = mrEM_MLE_G(sData,sKsi_tT,sNewG);
    sNewPhi = 1-sNewParam(2);
    sNewC = sNewParam(1);
    sNewSigma2 = sNewParam(3);
    sNewParam = [sNewPhi, sNewC, sNewSigma2, sNewG,  nan];

    % Base regime
    bNewG = fminsearch(@(g) mle_g(bData,bKsi_tT,g),oldParam(1,4),optimset('Display','off'));
    bNewParam = mrEM_MLE_G(bData,bKsi_tT,bNewG);
    bNewPhi = 1-bNewParam(2);
    bNewC = bNewParam(1);
    bNewSigma2 = bNewParam(3);
    bNewParam = [bNewPhi, bNewC, bNewSigma2, bNewG, nan];
    newParam = [bNewParam; sNewParam];
else
    % Use known g_1 and g_2
    % Spike regime
    sNewParam = mrEM_MLE_G(sData,sKsi_tT,oldParam(2,4));
    sNewPhi = 1-sNewParam(2);
    sNewC = sNewParam(1);
    sNewSigma2 = sNewParam(3);
    sNewParam = [sNewPhi, sNewC, sNewSigma2, oldParam(2,4), nan];

    % Base regime
    bNewParam = mrEM_MLE_G(bData,bKsi_tT,oldParam(1,4));
    bNewPhi = 1-bNewParam(2);
    bNewC = bNewParam(1);
    bNewSigma2 = bNewParam(3);
    bNewParam = [bNewPhi, bNewC, bNewSigma2, oldParam(1,4), nan];
    newParam = [bNewParam; sNewParam];
end

%==========================================================================
% mrs_EM_MLE auxiliary functions 
%==========================================================================

function LogL = mle_g(Data,KsiT_t,g)
beta = (sum(KsiT_t.*Data(1:end-1).*abs(Data(1:end-1)).^(-2*g).*( Data(2:end) - Data(1:end-1))) -...
    sum( KsiT_t ./ abs(Data(1:end-1)).^(2*g) .* ( Data(2:end) - Data(1:end-1) ) ) ./ sum ( KsiT_t ./ abs(Data(1:end-1)).^(2*g) ).*...
    sum(KsiT_t .*Data(1:end-1).*abs(Data(1:end-1)).^(-2*g) ))./...
    (sum(KsiT_t.*Data(1:end-1).*abs(Data(1:end-1)).^(-2*g)).* sum(KsiT_t.*Data(1:end-1).*abs(Data(1:end-1)).^(-2*g)) ./ sum( KsiT_t ./ abs(Data(1:end-1)).^(2*g) ) - sum(abs(Data(1:end-1)).^(2-2*g).*KsiT_t)  );
alpha = sum( KsiT_t ./ abs(Data(1:end-1)).^(2*g) .* ( Data(2:end) - Data(1:end-1) ))./sum(KsiT_t./abs(Data(1:end-1)).^(2*g))...
    + beta * sum(Data(1:end-1).*abs(Data(1:end-1)).^(-2*g).*KsiT_t )  ./ sum(KsiT_t./abs(Data(1:end-1)).^(2*g));
sigma2 = sum(KsiT_t./abs(Data(1:end-1)).^(2*g).*(Data(2:end)-Data(1:end-1)- alpha + beta*Data(1:end-1)).^2)./sum(KsiT_t);
LogL = abs(-sum(KsiT_t.*log(abs(Data(1:end-1))))+sum(log(abs(Data(1:end-1))).*(Data(2:end)-(1-beta).*Data(1:end-1)-alpha).^2./sigma2./abs(Data(1:end-1)).^(2*g).*KsiT_t));

function mrParam = mrEM_MLE_G(Data,KsiT_t,g)
beta = (sum(KsiT_t.*Data(1:end-1).*abs(Data(1:end-1)).^(-2*g).*( Data(2:end) - Data(1:end-1))) -...
    sum( KsiT_t ./ abs(Data(1:end-1)).^(2*g) .* ( Data(2:end) - Data(1:end-1) ) ) ./ sum ( KsiT_t ./ abs(Data(1:end-1)).^(2*g) ).*...
    sum(KsiT_t .*Data(1:end-1).*abs(Data(1:end-1)).^(-2*g) ))./...
    (sum(KsiT_t.*Data(1:end-1).*abs(Data(1:end-1)).^(-2*g)).* sum(KsiT_t.*Data(1:end-1).*abs(Data(1:end-1)).^(-2*g)) ./ sum( KsiT_t ./ abs(Data(1:end-1)).^(2*g) ) - sum(abs(Data(1:end-1)).^(2-2*g).*KsiT_t)  );
alpha = sum( KsiT_t ./ abs(Data(1:end-1)).^(2*g) .* ( Data(2:end) - Data(1:end-1) ))./sum(KsiT_t./abs(Data(1:end-1)).^(2*g))...
    + beta * sum(Data(1:end-1).*abs(Data(1:end-1)).^(-2*g).*KsiT_t )  ./ sum(KsiT_t./abs(Data(1:end-1)).^(2*g));
sigma2 = sum(KsiT_t./abs(Data(1:end-1)).^(2*g).*(Data(2:end)-Data(1:end-1)- alpha + beta*Data(1:end-1)).^2)./sum(KsiT_t);
mrParam = [alpha,beta,sigma2];

%==========================================================================
% mrs_EM_Smoother 
%==========================================================================

function [Ksi_tT,newP,newKsi_t1t_10,LogL] = mrs_EM_Smoother(Data,oldParam,oldP,oldKsi_t1t_10)

N = length(Data);
% Cond. dens. for the base regime (first column) and spike regime (second
% column)
sOldPhi = oldParam(2, 1);
sOldC = oldParam(2, 2);
sOldSigma2 = oldParam(2, 3);   
bOldPhi = oldParam(1, 1);
bOldC = oldParam(1, 2);
bOldSigma2 = oldParam(1, 3); 
Modelb = oldParam(1,4);
Model = oldParam(2,4);
Eta = [nan,nan;normpdf(Data(2:end), bOldPhi * Data(1:end-1) + bOldC, sqrt(bOldSigma2).*abs(Data(1:end-1)).^Modelb),...
      normpdf(Data(2:end), sOldPhi * Data(1:end-1) + sOldC, sqrt(sOldSigma2).*abs(Data(1:end-1)).^Model)];    
    
% Cond. prob. that t-th observation was generated by the base regime 
% (first column) or the spike regime (second column) based on the 
% knowledge of data up to time t
Ksi_tt = [oldKsi_t1t_10; zeros(N-1, 2)];

% Cond. prob. that (t+1)-st observation was generated by the base regime
% (first column) or the spike regime (second column) based on the 
% knowledge of data up to time t
Ksi_t1t = [oldKsi_t1t_10; zeros(N-1, 2)];

% Cond. prob. that t-th observation was generated by the base regime 
% (first column) or the spike regime (second column) based on the knowledge 
% of the whole data set (data up to time T), i.e. smoothed inferences
Ksi_tT = [oldKsi_t1t_10; zeros(N-1, 2)];

for row = 2:N
    Ksi_tt(row, :) = Eta(row, :) .* Ksi_t1t(row-1, :) / sum( Eta(row, :) .* Ksi_t1t(row-1, :) );
    Ksi_t1t(row, :) = Ksi_tt(row, :) * oldP;
end
Ksi_tT(end, :) = Ksi_tt(end, :);
for row = (N-1):-1:2
    Ksi_tT(row, :) = Ksi_tt(row, :) .* ( ( Ksi_tT(row+1, :) ./ Ksi_t1t(row, :) ) * oldP' );
end
    
% Log-likelihood
Aux = Eta(2:end,:) .* Ksi_tT(2:end,:);
LogL = sum(log(Aux(:,1) + Aux(:,2)));
  
% Transition probabilities
p_11 = sum( oldP(1,1) * Ksi_tT(3:end, 1) .* Ksi_tt(2:end-1, 1) ./ Ksi_t1t(2:end-1, 1) ) / sum( Ksi_tT(2:end-1, 1) );
p_22 = sum( oldP(2,2) * Ksi_tT(3:end, 2) .* Ksi_tt(2:end-1, 2) ./ Ksi_t1t(2:end-1, 2) ) / sum( Ksi_tT(2:end-1, 2) );
newP = [p_11, 1-p_11; 1-p_22, p_22];

newKsi_t1t_10 = Ksi_tT(2,:);

% Janczura, Weron: Matlab function to estimate parameters of a markov
% regime-switching (MRS) model with 2 independent regimes

function [Ksi_tT,Param,P,Ksi_t1t_10,LogL]=mrs2ir_est(Data,Model,g,qS,display,startParam,startP,startKsi_t1t_10)
%MRS2IR_EST Estimates parameters of a MRS model with 2 independent regimes. 
%   [KSI_TT,PARAM,P]=MRS2IR_EST(DATA,MODEL) returns smoothed inferences 
%   KSI_TT, estimated parameters PARAM and transition matrix P of a Markov 
%   regime-switching (MRS) model with 2 independent regimes: 
%    (i) an AR(1) process in the base regime
%        X(t+1)=phi_1*X(t)+c_1+sigma_1*|X(t)|^g_1*N(0,1),
%   (ii) and a Gaussian (MODEL='G'), shifted-lognormal ('LN'), 
%        shifted-Pareto ('P'), Weibull ('W') or exponential ('E') 
%        distributed spike regime,
%   fitted to time series DATA. 
%   The first column (KSI_TT) or row (PARAM, P) contains results for the
%   base regime and the second column/row for the spike regime.
%   [KSI_TT,PARAM,P,KSI_T1T_10,LOGL]=MRS2IR_EST(DATA,MODEL) additionally 
%   returns probabilities KSI_T1T_10 classifying the first observation to
%   one of the regimes and log-likelihood LOGL of the fitted model.
%
%   [...]=MRS2IR_EST(DATA,MODEL,G) allows to specify the known value of g_1
%   (default value G=[]). If G is not given (or G=[]) g_1 is estimated, 
%   otherwise a numerically less demanding algorithm is used for estimating 
%   the remaining parameters.
%
%   MRS2IR_EST(DATA,MODEL,G,QS) allows to specify a shift QS of the 
%   shifted-lognormal or shifted-Pareto distribution defining the spike 
%   regime (default value: QS=median(DATA)).
%
%   MRS2IR_EST(DATA,MODEL,G,QS,DISPLAY,STARTPARAM,STARTP,STARTKSIT1T_10) 
%   allows to specify initial (model-dependent) parameter estimates 
%   STARTPARAM (a 2x5 vector), initial transition matrix STARTP and initial 
%   estimates STARTKSIT1T_10 of probabilities classifying the first 
%   observation. Default values for the latter three parameters are:
%   	switch MODEL
%           case 'G'
%               STARTPARAM = [0.3, 15, 1, 0, NaN; NaN, 3, 0.5, NaN, NaN];
%           case 'LN'
%               STARTPARAM = [0.3, 15, 1, 0, NaN; NaN, 3, 1, qS, NaN];
%           case 'P'
%               STARTPARAM = [0.3, 15, 1, 0, NaN; NaN, NaN, NaN, 2, qS];
%           case 'W'
%               STARTPARAM = [0.3, 15, 1, 0, NaN; NaN, NaN, NaN, 1, 1];
%           case 'E'
%               STARTPARAM = [0.3, 15, 1, 0, NaN; NaN, NaN, NaN, NaN, 1];   
%       end
%       STARTP = [0.8, 0.2; 0.6, 0.4];
%       STARTKSIT1T_10 = [0.6, 0.4];
%   The fifth pameter, DISPLAY (default DISPLAY=1), is a flag which defines 
%   whether the calibration results are dispayed in the command window (1) 
%   or not (0).
%   
%   Examples:
%     1. Param = [0.2,10,0.01,1,nan;nan,20,2,0,nan]; P=[0.8,0.2;0.5,0.5];
%        [Y,S] = mrs2ir_sim(P,Param,'G',10,[1,0],2000,1);
%        [Ksi_tT,Param,P] = mrs2ir_est(Y,'G',1);
%        mrs2_plot(Ksi_tT,Y)
%
%     2. Param = [0.3,15,10,0,nan;nan,3,1,30,nan]; P=[0.8,0.2;0.5,0.5]; 
%        [Y,S] = mrs2ir_sim(P,Param,'LN',20,[1,0],2000,1);
%        [Ksi_tT,Param,P] = mrs2ir_est(Y,'LN',[],30); 
%        mrs2_plot(Ksi_tT,Y) 
% 
%   See also MRS2IR_SIM, MRS2_PLOT
% 

if nargin < 8
    startKsi_t1t_10 = [0.6, 0.4];
end
if nargin < 7
    startP = [0.8, 0.2; 0.6, 0.4];
end
if nargin<3
    g = [];
elseif ~isempty(g)
    % Assign known values of g_1
    startParam(1,4) = g(1);
end
if nargin<4
    % Spike distribution shift (lognormal, Pareto)
    qS = median(Data);
end
if nargin < 6
	switch Model
        case 'G'
            % Gaussian distribution
            startParam = [0.3, 15, 1, 0, NaN; NaN, 3, 0.5, NaN, NaN];
        case 'LN'
            % shifted-lognormal distribution
            startParam = [0.3, 15, 1, 0, NaN; NaN, 3, 1, qS, NaN];
        case 'P'
            % shifted-Pareto distribution
            startParam = [0.3, 15, 1, 0, NaN; NaN, NaN, NaN, 2, qS];
        case 'W'
            % Weibull distribution
            startParam = [0.3, 15, 1, 0, NaN; NaN, NaN, NaN, 1, 1];            
        case 'E'
            % Exponential distribution
            startParam = [0.3, 15, 1, 0, NaN; NaN, NaN, NaN, NaN, 1];  
    end
    T = median(Data);
    s = std(diff(Data));
    Pind = (Data>(T+s));
    startParam = mrs_EM_MLE(Data, Model, [1-Pind,Pind], [1-Pind,Pind], startParam, qS, g);    
end
if nargin<5
    display = 1;
end

% Initial estimation step 
[newKsi_tT, newKsi_tt, newP, newKsi_t1t_10, newLogL] = mrs_EM_Smoother(Data, Model, startParam, startP, startKsi_t1t_10, qS,g);   
newParam = mrs_EM_MLE(Data, Model, newKsi_tT,newKsi_tt, startParam, qS, g);

% Estimation loop 
oldLogL = newLogL;
iteration = 1;
dif = 1;
while and(dif > 10^(-8), iteration < 100)
    oldP = newP;
    oldKsi_t1t_10 = newKsi_t1t_10;
    oldParam = newParam;
    oldKsi_tT = newKsi_tT;
    oldLogL = newLogL;
    [newKsi_tT, newKsi_tt, newP, newKsi_t1t_10, newLogL] = mrs_EM_Smoother(Data, Model, oldParam, oldP, oldKsi_t1t_10, qS,g);    
    newParam = mrs_EM_MLE(Data, Model, newKsi_tT, newKsi_tt, oldParam, qS, g);
    
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
    P(ind)=0;
end
Ksi_t1t_10 = newKsi_t1t_10;
LogL = newLogL;
% Display results in the command window
if display
    Summary = mrs_Summary(Model, Param, P, LogL, qS);
    disp(['Gamma=' sprintf('%0.2g',Param(1,4)) ' for base regime']);
    switch Model
        case 'G'
            % Gaussian distribution
            disp([' ']);
            disp(['Two state regime switching model with Gaussian spikes, LogL=' num2str(LogL)]);
            disp([' ']);
            disp(['regime   phi_i   c_i      sigma^2_i  E(Y_{t,i})  Var(Y_{t,i})  q_ii     P(R=i)']);
            disp(['base    ' sprintf('%7.5f',Summary(1,1)) '  ' sprintf('%7.5f',Summary(1,2)) '  ' sprintf('%7.5f',Summary(1,3)) '    ' ...
            sprintf('%7.5f',Summary(1,4)) '     ' sprintf('%7.5f',Summary(1,5)) '       ' sprintf('%7.5f',Summary(1,6)) '  ' sprintf('%7.5f',Summary(1,7))]);
            disp(['spike   -        ' sprintf('%7.5f',Summary(2,2)) '  ' sprintf('%7.5f',Summary(2,3)) '    ' ...
            sprintf('%7.5f',Summary(2,4)) '     ' sprintf('%7.5f',Summary(2,5)) '       ' sprintf('%7.5f',Summary(2,6)) '  ' sprintf('%7.5f',Summary(2,7))]);    
        case 'LN'
            % log-normal distribution
            disp([' ']);
            disp(['Two state regime switching model with lognormal spikes, LogL=' num2str(LogL)]);        
            disp([' ']);
            disp(['regime   phi_i   c_i      sigma^2_i  E(Y_{t,i})  Var(Y_{t,i})  q_ii     P(R=i)']);
            disp(['base    ' sprintf('%7.5f',Summary(1,1)) '  ' sprintf('%7.5f',Summary(1,2)) '  ' sprintf('%7.5f',Summary(1,3)) '    ' ...
                sprintf('%7.5f',Summary(1,4)) '     ' sprintf('%7.5f',Summary(1,5)) '       ' sprintf('%7.5f',Summary(1,6)) '  ' sprintf('%7.5f',Summary(1,7))]);
            disp(['spike           '  sprintf('%7.5f',Summary(2,2)) '  ' sprintf('%7.5f',Summary(2,3)) '    ' ...
                sprintf('%7.5f',Summary(2,4)+qS) '     ' sprintf('%7.5f',Summary(2,5)) '       ' sprintf('%7.5f',Summary(2,6)) '  ' sprintf('%7.5f',Summary(2,7))]);    
        case 'P'
            % Pareto distribution
            disp([' ']);        
            disp(['Two state regime switching model with Pareto spikes, LogL=' num2str(LogL)]);                
            disp([' ']);        
            disp(['regime   phi_i   c_i / a  s^2_i / k  E(Y_{t,i})  Var(Y_{t,i})  q_ii     P(R=i)']);
            disp(['base    ' sprintf('%7.5f',Summary(1,1)) '  ' sprintf('%7.5f',Summary(1,2)) '  ' sprintf('%7.5f',Summary(1,3)) '    ' ...
                sprintf('%7.5f',Summary(1,4)) '     ' sprintf('%7.5f',Summary(1,5)) '       ' sprintf('%7.5f',Summary(1,6)) '  ' sprintf('%7.5f',Summary(1,7))]);
            disp(['spike   -        ' sprintf('%7.5f',Summary(2,2)) '  ' sprintf('%7.5f',Summary(2,3)) '    ' ...
                sprintf('%7.5f',Summary(2,4)) '     ' sprintf('%7.5f',Summary(2,5)) '       ' sprintf('%7.5f',Summary(2,6)) '  ' sprintf('%7.5f',Summary(2,7))]);    
        case 'W'
            % Weibull distribution
            disp([' ']);        
            disp(['Two state regime switching model with Weibull spikes, LogL=' num2str(LogL)]);                
            disp([' ']);        
            disp(['regime   phi_i   c_i / a  s^2_i / b  E(Y_{t,i})  Var(Y_{t,i})  q_ii     P(R=i)']);
            disp(['base    ' sprintf('%7.5f',Summary(1,1)) '  ' sprintf('%7.5f',Summary(1,2)) '  ' sprintf('%7.5f',Summary(1,3)) '    ' ...
                sprintf('%7.5f',Summary(1,4)) '     ' sprintf('%7.5f',Summary(1,5)) '       ' sprintf('%7.5f',Summary(1,6)) '  ' sprintf('%7.5f',Summary(1,7))]);
            disp(['spike   -        ' sprintf('%7.5f',Summary(2,2)) '  ' sprintf('%7.5f',Summary(2,3)) '    ' ...
                sprintf('%7.5f',Summary(2,4)) '     ' sprintf('%7.5f',Summary(2,5)) '       ' sprintf('%7.5f',Summary(2,6)) '  ' sprintf('%7.5f',Summary(2,7))]);    
        case 'E'
            % Exponential distribution
            disp([' ']);        
            disp(['Two state regime switching model with Exponential spikes, LogL=' num2str(LogL)]);                
            disp([' ']);        
            disp(['regime   phi_i   c_i / a  s^2_i / l  E(Y_{t,i})  Var(Y_{t,i})  q_ii     P(R=i)']);
            disp(['base    ' sprintf('%7.5f',Summary(1,1)) '  ' sprintf('%7.5f',Summary(1,2)) '  ' sprintf('%7.5f',Summary(1,3)) '    ' ...
                sprintf('%7.5f',Summary(1,4)) '     ' sprintf('%7.5f',Summary(1,5)) '       ' sprintf('%7.5f',Summary(1,6)) '  ' sprintf('%7.5f',Summary(1,7))]);
            disp(['spike   -        ' '-      ' '  ' sprintf('%7.5f',Summary(2,3)) '    ' ...
                sprintf('%7.5f',Summary(2,4)) '     ' sprintf('%7.5f',Summary(2,5)) '       ' sprintf('%7.5f',Summary(2,6)) '  ' sprintf('%7.5f',Summary(2,7))]);    
    end
end

%=========================================================================
% Internally used routine(s)
%=========================================================================

%=========================================================================
% mrs_Summary
%=========================================================================

function Summary = mrs_Summary(Model, Param, P, logL, qS)

% Mean and variance for the base regime
bParam = Param(1,:);
bPhi = bParam(1);
bC = bParam(2);
bSigma2 = bParam(3);
bMean = bC / (1-bPhi);
if bParam(4)==0
    bVar = bSigma2 / 2/(1-bPhi);
else
    bVar = NaN;
end

% Mean and variance for the spike regime
sParam = Param(2,:);
[sMean, sVar] = mrs_Moments(Model, sParam);

% Unconditional probabilities
p_11 = P(1,1);
p_22 = P(2,2);
P_Rt_1 = (1-p_22) / (2-p_11-p_22);
P_Rt_2 = (1-p_11) / (2-p_11-p_22);

% Summary matrix
switch Model
    case 'G'
        %Gaussian distribution
        sC = sParam(2);
        sSigma2 = sParam(3);
        Summary = [ bPhi,   bC,     bSigma2,   bMean,  bVar,   p_11,   P_Rt_1;
                    nan,    sC,     sSigma2,   sMean,  sVar,   p_22,   P_Rt_2;
                    logL,   nan,    nan,       nan,    nan,    nan,    nan];
    case 'LN'
        sC = sParam(2);
        sSigma2 = sParam(3);
        Summary = [ bPhi,   bC,     bSigma2,   bMean,  bVar,   p_11,   P_Rt_1;
                    nan,    sC,     sSigma2,   sMean,  sVar,   p_22,   P_Rt_2;
                    logL,   nan,    nan,       nan,    nan,    nan,    nan];
    case 'P'
        %Pareto distribution
        sA = sParam(4);
        sK = sParam(5);
        Summary = [ bPhi,   bC,     bSigma2,   bMean,  bVar,   p_11,   P_Rt_1;
                    nan,    sA,     sK,        sMean,  sVar,   p_22,   P_Rt_2;
                    logL,   nan,    nan,       nan,    nan,    nan,    nan];
    case 'W'
        %Weibull distribution
        sA = sParam(4);
        sB = sParam(5);
        Summary = [ bPhi,   bC,     bSigma2,   bMean,  bVar,   p_11,   P_Rt_1;
                    nan,    sA,     sB,        sMean,  sVar,   p_22,   P_Rt_2;
                    logL,   nan,    nan,       nan,    nan,    nan,    nan];    
    case 'E'
        %Exponential distribution
        sB = sParam(5);
        Summary = [ 1-bPhi,   bC,     bSigma2,   bMean,  bVar,   p_11,   P_Rt_1;
                    nan,      nan,    sB,        sMean,  sVar,   p_22,   P_Rt_2;
                    logL,     nan,    nan,       nan,    nan,    nan,    nan];
end

%=========================================================================
% mrs_Moments
%=========================================================================

function [sMean, sVar] = mrs_Moments(Model, sParam)

if strcmp(Model, 'G') == 1
    %Gaussian distribution:
    sC = sParam(2);
    sSigma2 = sParam(3);
    sMean = sC;
    sVar = sSigma2;
elseif strcmp(Model, 'LN') == 1
    %log-normal distribution
    sC = sParam(2);
    sSigma2 = sParam(3);
    sMean = exp(sC + sSigma2/2);
    sVar = exp(2*sC+sSigma2)*(exp(sSigma2)-1); 
elseif strcmp(Model, 'P') == 1
    %Pareto distribution
    sA = sParam(4);
    sK = sParam(5);   
    sMean = (sA*sK)/(sA-1);
    sVar = (sA*sK*sK)/((sA-1)*(sA-1)*(sA-2));
    if sA<=2
        sVar = Inf;
    end
    if sA<=1
        sMean = Inf;
    end    
elseif strcmp(Model, 'W') == 1
    %Weibull distribution
    sA = sParam(4);
    sB = sParam(5);
    sMean = sA^(-1/sB)*gamma(1+1/sB);
    sVar = sA^(-2/sB)*(gamma(1+2/sB)-(gamma(1+1/sB))^2);
    if a<2
        sVar = nan;
    end
    if a<1
        sMean = nan;
    end
elseif strcmp(Model, 'E') == 1
    %Exponential distribution
    sB = sParam(5);
    sMean = sB;
    sVar = sB^2;
end

%=========================================================================
% mrs_EM_MLE
%=========================================================================

function newParam = mrs_EM_MLE(Data, Model, Ksi_tT, Ksi_tt, oldParam, qS, g)
% MRS_EM_MLE ML estimation of regime processes parameters. Note that 
% Pareto distribution uses moment estimators.

sData = Data(2:end);
sKsi_tT = Ksi_tT(2:end, 2);
bKsi_tT = Ksi_tT(2:end, 1);
% Spike regime 
switch Model
    case 'G'
        %Gaussian distribution
        sNewC = sum((sData).*sKsi_tT) / sum(sKsi_tT);
        sNewSigma2 = sum(((sData)-sNewC).^2 .* sKsi_tT) / sum(sKsi_tT);
        sNewParam = [nan, sNewC, sNewSigma2, nan, nan];
    case 'LN'
        %log-normal distribution
        ind = find(sData>qS);
        sNewC = sum(log(sData(ind)-qS).*sKsi_tT(ind)) / sum(sKsi_tT(ind));
        sNewSigma2 = sum((log(sData(ind)-qS)-sNewC).^2 .* sKsi_tT(ind)) /sum(sKsi_tT(ind));
        sNewParam = [nan,sNewC, sNewSigma2, qS, nan];
    case 'P'
        %Pareto distribution
        sMin=qS;
        sNewK = sMin;
        sNewA = sum(sKsi_tT)/sum(sKsi_tT .* (log(sData)-log(sNewK)));
        sNewParam = [nan, nan, nan, sNewA, sNewK];       
    case 'W'
        %Weibull distribution
        sOldA = oldParam(2,4);
        sOldB = oldParam(2,5);
        sNewParam = fsolve(@weibullEM_MLE, [sOldA, sOldB], optimset('Display', 'off'), sData, sKsi_tT);
        sNewA = sNewParam(1);
        sNewB = sNewParam(2);
        sNewParam = [nan, nan, nan, sNewA, sNewB];
    case 'E'
        %Exponential distribution
        sNewParam =sum(sData.*sKsi_tT) / sum(sKsi_tT);
        sNewB = sNewParam;
        sNewParam = [nan, nan, nan, nan, sNewB];
end

bData = Data(1:end);
EData=Data(1:end);

% Base regime ML estimation
N = length(Data);
bOldPhi = oldParam(1,1);
bOldC = oldParam(1,2);
for row = 2:N       
    if Ksi_tT(row, 1)<1
        EData(row) = bData(row).*Ksi_tt(row,1) + (bOldPhi.*EData(row-1)+bOldC).*Ksi_tt(row,2);
    end    
end
if isempty(g)
    % Gamma estimated from data
    bNewG = fminsearch(@(g) mle_g(bData,EData,bKsi_tT,g),oldParam(1,4),optimset('Display','off'));
    bNewParam = mrEM_MLE_G(bData,EData,bKsi_tT,bNewG);
    bNewPhi = 1-bNewParam(2);
    bNewC = bNewParam(1);
    bNewSigma2 = bNewParam(3);
    bNewParam = [bNewPhi,bNewC,bNewSigma2,bNewG,nan];
    newParam = [bNewParam; sNewParam];
else    
    % Gamma known
    bNewParam = mrEM_MLE_G(bData,EData,bKsi_tT,g);
    bNewPhi = 1-bNewParam(2);
    bNewC = bNewParam(1);
    bNewSigma2 = bNewParam(3);
    bNewParam = [bNewPhi,bNewC,bNewSigma2,g,nan];
    newParam = [bNewParam; sNewParam];
end

%==========================================================================
% mrs_EM_MLE auxiliary functions 
%==========================================================================

function LogL = mle_g(Data,EData,KsiT_t,g)
beta = (sum(KsiT_t.*EData(1:end-1).*abs(EData(1:end-1)).^(-2*g).*( Data(2:end) - EData(1:end-1))) -...
    sum( KsiT_t ./ abs(EData(1:end-1)).^(2*g) .* ( Data(2:end) - EData(1:end-1) ) ) ./ sum ( KsiT_t ./ abs(EData(1:end-1)).^(2*g) ).*...
    sum(KsiT_t .*EData(1:end-1).*abs(EData(1:end-1)).^(-2*g) ))./...
    (sum(KsiT_t.*EData(1:end-1).*abs(EData(1:end-1)).^(-2*g)).* sum(KsiT_t.*EData(1:end-1).*abs(EData(1:end-1)).^(-2*g)) ./ sum( KsiT_t ./ abs(EData(1:end-1)).^(2*g) ) - sum(abs(EData(1:end-1)).^(2-2*g).*KsiT_t)  );
alpha = sum( KsiT_t ./ abs(EData(1:end-1)).^(2*g) .* ( Data(2:end) - EData(1:end-1) ))./sum(KsiT_t./abs(EData(1:end-1)).^(2*g))...
    + beta * sum(EData(1:end-1).*abs(EData(1:end-1)).^(-2*g).*KsiT_t )  ./ sum(KsiT_t./abs(EData(1:end-1)).^(2*g));
sigma2 = sum(KsiT_t./abs(EData(1:end-1)).^(2*g).*(Data(2:end)-EData(1:end-1)- alpha + beta*EData(1:end-1)).^2)./sum(KsiT_t);
LogL = abs(-sum(KsiT_t.*log(abs(EData(1:end-1))))+sum(log(abs(EData(1:end-1))).*(Data(2:end)-(1-beta).*EData(1:end-1)-alpha).^2./sigma2./abs(EData(1:end-1)).^(2*g).*KsiT_t));

function mrParam = mrEM_MLE_G(Data,EData,KsiT_t,g)
beta = (sum(KsiT_t.*EData(1:end-1).*abs(EData(1:end-1)).^(-2*g).*( Data(2:end) - EData(1:end-1))) -...
    sum( KsiT_t ./ abs(EData(1:end-1)).^(2*g) .* ( Data(2:end) - EData(1:end-1) ) ) ./ sum ( KsiT_t ./ abs(EData(1:end-1)).^(2*g) ).*...
    sum(KsiT_t .*EData(1:end-1).*abs(EData(1:end-1)).^(-2*g) ))./...
    (sum(KsiT_t.*EData(1:end-1).*abs(EData(1:end-1)).^(-2*g)).* sum(KsiT_t.*EData(1:end-1).*abs(EData(1:end-1)).^(-2*g)) ./ sum( KsiT_t ./ abs(EData(1:end-1)).^(2*g) ) - sum(abs(EData(1:end-1)).^(2-2*g).*KsiT_t)  );
alpha = sum( KsiT_t ./ abs(EData(1:end-1)).^(2*g) .* ( Data(2:end) - EData(1:end-1) ))./sum(KsiT_t./abs(EData(1:end-1)).^(2*g))...
    + beta * sum(EData(1:end-1).*abs(EData(1:end-1)).^(-2*g).*KsiT_t )  ./ sum(KsiT_t./abs(EData(1:end-1)).^(2*g));
sigma2 = sum(KsiT_t./abs(EData(1:end-1)).^(2*g).*(Data(2:end)-EData(1:end-1)- alpha + beta*EData(1:end-1)).^2)./sum(KsiT_t);
mrParam = [alpha,beta,sigma2];

function F = weibullEM_MLE(wParam, wData, wKsi_tT)
% Weibull distribution likelihood function
wA = wParam(1);
wB = wParam(2);
F = [sum((1/wA - wData.^wB).*wKsi_tT);
    sum((1/wB + log(wData) - wA * log(wData) .* wData.^wB).*wKsi_tT)];

%==========================================================================
% mrs_EM_Smoother 
%==========================================================================

function [Ksi_tT, Ksi_tt, newP, newKsi_t1t_10, LogL] = mrs_EM_Smoother(Data, Model, oldParam, oldP, oldKsi_t1t_10, qS,g)

N = length(Data);
% Old base regime parameters
bOldPhi = oldParam(1, 1);
bOldC = oldParam(1, 2);
bOldSigma2 = oldParam(1, 3);

% Cond. dens. for the spike regime 
switch Model
    case 'G'
        %Gaussian distribution
        sOldC = oldParam(2,2);
        sOldSigma2 = oldParam(2,3);
        EtaS = [nan; normpdf(Data(2:end), sOldC, sqrt(sOldSigma2))];   
    case 'LN'
        sOldC = oldParam(2,2);
        sOldSigma2 = oldParam(2,3);
        EtaS = [nan; lognSpdf(Data(2:end), sOldC, sqrt(sOldSigma2),qS)];
    case 'P'
        %Pareto distribution
        sOldA = oldParam(2,4);
        sOldK = oldParam(2,5);
        EtaS = [nan;  paretopdf(Data(2:end), sOldA, sOldK)];
    case 'W'
        %Weibull distribution
        sOldA = oldParam(2,4);
        sOldB = oldParam(2,5);
        EtaS = [nan;  weibullpdf(Data(2:end), sOldA, sOldB)];      
    case 'E'
        %Exponential distribution
        sOldB = oldParam(2,5);
        EtaS = [nan;  exppdf(Data(2:end), sOldB)];     
end
    
% Cond. dens. for the base regime 
% Cond. prob. that t-th observation was generated by the mean-reverting 
% regime (first column) or the spike regime (second column) based on the 
% knowledge of data up to time t
Ksi_tt = [oldKsi_t1t_10; zeros(N-1, 2)];

% Cond. prob. that (t+1)-st observation was generated by the mean-reverting 
% regime (first column) or the spike regime (second column) based on the 
% knowledge of data up to time t
Ksi_t1t = [oldKsi_t1t_10; zeros(N-1, 2)];

% Cond. prob. that t-th observation was generated by the mean-reverting regime 
% (first column) or the spike regime (second column) based on the knowledge 
% of the whole data set (data up to time T), i.e. smoothed inferences
Ksi_tT = [oldKsi_t1t_10; zeros(N-1, 2)];

% Predictions about unobserved data from the base regime
EData = Data(1:end);
% Conditional density for the base regime
EtaB = zeros(N,1);
if isempty(g)
    % Gamma estimated from data
    Modelb = oldParam(1,4);
else
    % Gamma known
    Modelb = g;
end
for row = 2:N       
    EtaB(row) = normpdf(Data(row), bOldPhi .*EData(row-1) + bOldC, sqrt(bOldSigma2).*abs(EData(row-1)).^Modelb);
    Eta = [EtaB(row),EtaS(row)];
    Ksi_tt(row, :) = Eta.* Ksi_t1t(row-1, :) / sum( Eta .* Ksi_t1t(row-1, :) );
    Ksi_t1t(row, :) = Ksi_tt(row, :) * oldP;
    if Ksi_tt(row, 1)<1
        EData(row) = Data(row).*Ksi_tt(row,1)+(bOldPhi.*EData(row-1)+bOldC).*Ksi_tt(row,2);
    end    
end
Ksi_tT(end, :) = Ksi_tt(end, :);

for row = (N-1):-1:2
    Ksi_tT(row, :) = Ksi_tt(row, :) .* ( ( Ksi_tT(row+1, :) ./ Ksi_t1t(row, :) ) * oldP' );
end
%Log-likelihood
Eta = [EtaB,EtaS];
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

function Y = paretopdf(X, A, K)
n = length(X);
Y = zeros(n,1);
ind = find(X>K);
Y(ind) = (A*(K^A))./(X(ind).^(A+1));

function Y = weibullpdf(X, A, B)
if (A>0 & B>0 & sum(X<=0)==0)
    Y = A*B*(X.^(B-1)).*exp(-A*(X.^B));
else
    Y = nan;
end

function Y=lognSpdf(X, mu, sigma,c)
Y = zeros(length(X),1);
ind = find(X>c);
Y(ind) = lognpdf_mod(X(ind)-c,mu,sigma);

function Y=lognpdf_mod(X, mu, sigma)
Y = lognpdf(X, mu, sigma);
for i=1:length(X)
    if isnan(Y(i))
        Y(i) = 0;
    end
end

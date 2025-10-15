% Janczura, Weron: Matlab function to simulate trajectories of  a markov
% regime-switching (MRS) model with 2 independent regimes

function [Y,S]=mrs2ir_sim(P,Param,Model,StartValue,Ksi_0,T,N)
%MRS2IR_SIM Simulates trajectories of a MRS model with 2 independent regimes. 
%   MRS2IR_SIM(P,PARAM,MODEL,STARTVALUE,KSI_0,T,N) generates N trajectories 
%   of a Markov regime-switching (MRS) model with 2 independent regimes: 
%    (i) an AR(1) process in the base regime
%        X(t+1)=phi_1*X(t)+c_1+sigma_1*|X(t)|^g_1*N(0,1),
%   (ii) and a Gaussian (MODEL='G'), lognormal ('LN'), Pareto ('P'), 
%        Weibull ('W') or exponential ('E') distributed spike regime. 
%   P, PARAM and KSI_0 are the transition matrix, model parameters and 
%   probabilities classifying the first observation to one of the regimes, 
%   respectively, returned by MRS2IR_EST. STARTVALUE is the starting point 
%   for the simulated trajectories and T is the length of the trajectories.
%   [Y,S]=MRS2IR_SIM(P,PARAM,MODEL,STARTVALUE,T,N) returns also matrix S 
%   with ones, if the respective observation was generated from the spike 
%   regime, and zeros otherwise.
%   
%   Examples:
%     1. Param = [0.2,10,0.01,1,nan;nan,20,2,0,nan]; P=[0.8,0.2;0.5,0.5];
%        [Y,S] = mrs2ir_sim(P,Param,'G',10,[1,0],2000,1);
%
%     2. Param = [0.3,15,10,0,nan;nan,3,1,30,nan]; P=[0.8,0.2;0.5,0.5]; 
%        [Y,S] = mrs2ir_sim(P,Param,'LN',20,[1,0],2000,1);
%
%   See also MRS2IR_EST, MRS2_PLOT
%
%   Reference(s):
%   [1] J.Janczura, R.Weron (2011) Efficient estimation of Markov 
%   regime-switching models: An application to electricity spot prices. 
%   Working paper version available at: 
%   http://ideas.repec.org/p/wuu/wpaper/hsc1102.html   
%   [2] R.Weron (2007) Modeling and Forecasting Electricity Loads and 
%   Prices: A Statistical Approach, Wiley, Chichester. 

%   Written by Joanna Janczura (2010.08.17)
%   Revised by Rafal Weron (2011.10.03)

Y = zeros(T+1,N);
Y(1,:) = StartValue;
S = zeros(T+1,N);
ParamB = Param(1,:);
ParamS = Param(2,:);
YtB = g_sim(ParamB(1),ParamB(2),sqrt(ParamB(3)),ParamB(4),Y(1),T,1,N);
for j = 1:N
    r = rand<Ksi_0(2)+1;    
    if r==2
       S(1,j) = 1;
    end
    for i = 2:T+1
        if rand<P(r,1)
            Y(i,j) = YtB(i,j);
            r = 1;
        else
            Y(i,j) = YtS(ParamS,Model);
            r = 2;
            S(i,j) = 1;
        end     
    end
end

%=========================================================================
% Internally used routine(s)
%=========================================================================

function X = g_sim(phi,alpha,sigma,g,x0,T,dt,N)
% X=G_SIM(PHI,ALPHA,SIGMA,G,X0,T,DT,N) returns a matrix of N realizations
% of a mean reverting process defined by the following SDE: 
%   dX = (ALPHA - (1-PHI)*X) * DT + SIGMA * |X|^G * dB_t
% for t=0,DT,2*DT,...,T with an initial value of the process X0.

B = normrnd(0,1,T/dt,N);
X = zeros(T/dt+1,N);
X(1,:) = x0;
for i = 1:(T/dt)
   X(i+1,:) = X(i,:) + (alpha - (1-phi)*X(i,:))*dt + sigma*sqrt(dt).*abs(X(i,:)).^g.*B(i,:);
end

function Z = YtS(ParamS,Model)
switch Model
    case 'G'
        Z = normrnd(ParamS(2),sqrt(ParamS(3)));
    case 'LN'
        Z = lognrnd(ParamS(2),sqrt(ParamS(3)))+ParamS(4);
    case 'P'
        Z = ParamS(5)*(rand.^(-1/ParamS(4)));
    case 'W'
        Z = wblrnd(ParamS(4),ParamS(5));
    case 'E'
        Z = exprnd(ParamS(5));
end

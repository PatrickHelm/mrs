% Janczura, Weron: Matlab function ti simulate trajectories of a 2-regime
% parameter switching (PS) model

function [Y,S]=ps2r_sim(P,Param,StartValue,Ksi_0,T,N)
%PS2R_SIM Simulates trajectories of a 2-regime parameter switching (PS) model. 
%   PS2R_SIM(P,PARAM,STARTVALUE,KSI_0,T,N) generates N trajectories of 
%   a 2-regime parameter switching (PS) model, i.e. a 2-regime Markov 
%   regime-switching (MRS) model with both regimes driven by AR(1) 
%   processes of the form: X(t+1)=phi_i*X(t)+c_i+sigma_i*|X(t)|^g_i*N(0,1). 
%   P, PARAM and KSI_0 are the transition matrix, model parameters and 
%   probabilities classifying the first observation to one of the regimes, 
%   respectively, returned by PS2R_EST. STARTVALUE is the starting point 
%   for the simulated trajectories and T is the length of the trajectories.
%   [Y,S]=PS2R_SIM(P,PARAM,STARTVALUE,T,N) returns also matrix S with ones, 
%   if the respective observation was generated from the spike regime, and 
%   zeros otherwise.
%   
%   Example:
%       Param = [0.2,2,1,0;0.4,3,1,1]; P=[0.5,0.5;0.4,0.6];
%       [Y,S] = ps2r_sim(P,Param,3.5,[1,0],1000,1);
%
%   See also PS2R_EST
%
%   Reference(s):
%   [1] J.Janczura, R.Weron (2011) Efficient estimation of Markov 
%   regime-switching models: An application to electricity spot prices. 
%   Working paper version available at: 
%   http://ideas.repec.org/p/wuu/wpaper/hsc1102.html

%   Written by Joanna Janczura (2010.05.31)
%   Revised by Rafal Weron (2011.10.03)  

ParamB = Param(1,:);
ParamS = Param(2,:);

Y = zeros(T+1,N);
Y(1,:) = StartValue;
S = zeros(T+1,N);
for j = 1:N
    r = rand < Ksi_0(2)+1;    
    S(1,j) = (r==2);
    for i = 2:T+1
        if rand < P(r,1)
            Y(i,j) = ParamB(1)*Y(i-1,j) + ParamB(2) + ...
                normrnd(0,1,1,1).*sqrt(ParamB(3)).*abs(Y(i-1,j)).^ParamB(4);
            r = 1;
        else
            Y(i,j) = ParamS(1)*Y(i-1,j) + ParamS(2) + ...
                normrnd(0,1,1,1).*sqrt(ParamS(3)).*abs(Y(i-1,j)).^ParamS(4);
            r = 2;
            S(i,j) = 1;
        end     
    end
end
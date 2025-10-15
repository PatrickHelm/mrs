% Example Script MS_Regress_Fit.m based on MS_Regress toolbox of Perlin

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');
logRet=importdata('Example_Oil_10w.txt');  % load some Data.

dep=logRet(:,1);                    % Defining dependent variable from .mat file
constVec=ones(length(dep),1);       % Defining a constant vector in mean equation (just an example of how to do it)
indep=[constVec logRet(:,2:3)];     % Defining some explanatory variables
k=3;                                % Number of States
S=[1 0 0 0];                        % first: mean, last: variance
advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details

[Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model

rmpath('m_Files');
rmpath('data_Files'); 
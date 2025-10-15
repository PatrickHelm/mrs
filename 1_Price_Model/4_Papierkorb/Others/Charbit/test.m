clear all
addpath ../toolbox
% sp_three obtained with the program "xCEPS\excc.m" 
% 10 occurences of the number "3" in english
load ../xCEPS/sp_three 
cc_training=cc(:,:,2:end);
[Mu, Sigma, A, Pi, LL_training, ct] = hmmtrain(cc_training,Tr,4);


cc_test=cc(:,:,1);
LL_reco = hmmrecog(cc_test, Tr(end), Mu, Sigma, A, Pi);
LL_reco
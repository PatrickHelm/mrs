function MarkovModel = getMarkovModel(pDist)
%GETMARKOVMODEL Creates the transition probability matrix for R1 and R2 
%  MarkoModel=GETMARKOVMODEL(pDist) creates a 1-by-2 cell containg the 
%  transition probability matrix for R1 and R2 according to Gavirneni(2004).
%  PDist itself is a cell 1-by-2 cell containing the name of the model for
%  exchange rate fluctuations. This can either be
%  - 'random walk'
%  - 'mean-reverting'
%  - 'momentum'

MarkovModel=cell(1,2);
for R=1:2
    if strcmp(pDist{R}, 'RW')
        MarkovModel{R} = [0.5 0.5 0 0 0; 0.5 0 0.5 0 0 ; 0 0.5 0 0.5 0; 0 0 0.5 0 0.5; 0 0 0 0.5 0.5];
    elseif strcmp(pDist{R}, 'MR')
        MarkovModel{R} = [0.2 0.8 0 0 0; 0.1 0.3 0.6 0 0 ; 0 0.1 0.8 0.1 0; 0 0 0.6 0.3 0.1; 0 0 0 0.8 0.2];
	elseif strcmp(pDist{R}, 'MO')
        MarkovModel{R} = [0.8 0.2 0 0 0; 0.6 0.3 0.1 0 0 ; 0 0.5 0 0.5 0; 0 0 0.1 0.3 0.6; 0 0 0 0.2 0.8];
    end
end
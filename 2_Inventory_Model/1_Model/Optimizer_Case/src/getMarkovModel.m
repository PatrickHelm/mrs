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
        MarkovModel{R} = [0.83 0.17 0     0     0     0     0     0; 
                          0.13 0.73 0.13  0     0     0     0     0; 
                          0    0.15 0.55  0.15  0.15  0     0     0; 
                          0    0    0.17  0.40  0.43  0     0     0;
                          0    0    0     0.14  0.43  0.29  0.14  0;
                          0    0    0     0     0.38  0.50  0.00  0.13;
                          0    0    0     0     0     0.50  0.50  0;
                          0    0    0     0     0     0.33  0     0.67];
    elseif strcmp(pDist{R}, 'MR')
        MarkovModel{R} = [0    0    1     0     0     0     0     0; 
                          0    0.80 0.20  0     0     0     0     0; 
                          0    0.02 0.88  0.11  0     0     0     0; 
                          0    0    0.5   0.5   0     0     0     0;
                          0    0    1     0     0     0     0     0; 
                          0    0    1     0     0     0     0     0; 
                          0    0    1     0     0     0     0     0; 
                          0    0    1     0     0     0     0     0];
    end
end
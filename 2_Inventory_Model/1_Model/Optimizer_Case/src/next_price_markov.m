%nextPrice calculates the probabilities for all possible prices in t+1

function nP_markov = next_price_markov(markovmodel,prices, t, PI_markov)

NoP=numel(prices); 
nP_markov=zeros(NoP^t, NoP);
R1=markovmodel{1}; %probability transition matrix for R1
R2=markovmodel{2}; %probability transition matrix for R2

if t==1
    for p_t=1:NoP %p_t is price observation in period t
        for p_future=1:NoP %p_future = p_t+1
            nP_markov(p_t, p_future)=PI_markov(1)*R1(p_t, p_future)+PI_markov(2)*R2(p_t, p_future);
        end
    end
else
    PI_R1=PI_markov{t}{1}; %next period regime belief at time t for R1
    PI_R2=PI_markov{t}{2}; %next period regime belief at time t for R2
    j=1; %index for resulting matrix
    i=1; %index PI_R1 and PI_R2;
    while j<=NoP^t
        for p_t=1:NoP
           for p_future=1:NoP
              nP_markov(j, p_future)=PI_R1(i,p_t)*R1(p_t, p_future)+PI_R2(i,p_t)*R2(p_t, p_future);
              %nP_Markov(j, p_future)=PI_R1(p_t,i)*R1(i, p_future)+PI_R2(p_t,i)*R2(i, p_future);
           end
           j=j+1;
        end
        if mod(i,NoP)==0
            i=0;
        end
        i=i+1;
    end
end
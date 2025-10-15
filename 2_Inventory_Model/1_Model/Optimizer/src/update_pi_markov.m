%update_pi: Calculation of PI_markov

function PI_markov = update_pi_markov(markovmodel,prices,t,PI,k,solution_method)

NoP=numel(prices); %Number of Prices
PI_markov=cell(1,2);
R1=markovmodel{1}; %probability transition matrix for R1
R2=markovmodel{2}; %probability transition matrix for R2
if t==1 %for t=1 the belief update values are given
    PI_markov{1}=PI(1);
    PI_markov{2}=PI(2);
    return
end

PI_markov1=zeros(NoP^(t-1), NoP);
PI_markov2=zeros(NoP^(t-1), NoP);
PI_past_R1=PI{t-1}{1}; %belief updates in t-1 for R1
PI_past_R2=PI{t-1}{2}; %belief updates in t-1 for R2
if t==2 % In t=1, pi=[0.5 0.5]. To use only one loop for all t>=2, we have to blow up the size of the belief in period 1 by copying its value 5 (NoP) times
    PI_past_R1=repmat(PI_past_R1,1,NoP);
    PI_past_R2=repmat(PI_past_R2,1,NoP);
end
i=1;
j=1;
while j<=NoP^(t-1)
    for last_price=1:NoP %p_past is the price in the last period p_(t-1)
       for p_t=1:NoP %p_t is the price in the actual period t
           if strcmp(solution_method, 'CEC')  
                    PI_markov1(j,p_t)=PI{1}{1};
                    PI_markov2(j, p_t)=PI{1}{2};
           elseif strcmp(solution_method, 'R1')
                     PI_markov1(j,p_t)=1;
                     PI_markov2(j,p_t)=0;
           elseif strcmp(solution_method, 'R2')
                     PI_markov1(j,p_t)=0;
                     PI_markov2(j,p_t)=1;
           else 
          PI_markov1(j,p_t) = (PI_past_R1(i,last_price)*k(1,1)*R1(last_price, p_t)...
              +PI_past_R2(i,last_price)*k(1,2)*R2(last_price,p_t))... %enumerator for R1
              /(PI_past_R1(i,last_price)*R1(last_price, p_t)+PI_past_R2(i,last_price)*R2(last_price, p_t)); %denominator for R1
          PI_markov2(j, p_t) = (PI_past_R1(i,last_price)*k(2,1)*R1(last_price, p_t)... 
              +PI_past_R2(i,last_price)*k(2,2)*R2(last_price, p_t))...%enumerator for R2
              /(PI_past_R1(i,last_price)*R1(last_price, p_t)+(PI_past_R2(i,last_price)*R2(last_price, p_t))); %denominator for R1
          if isnan(PI_markov1(j, p_t)) %Division by zero -> NO belief update
              PI_markov1(j, p_t)=PI_past_R1(i,last_price);
          end
          if isnan(PI_markov2(j, p_t)) %Division by zero -> NO belief update
              PI_markov2(j, p_t)=PI_past_R2(i,last_price);
          end
           end %neu
       end
       j=j+1;
    end
    i=i+1;
end

PI_markov={PI_markov1, PI_markov2};


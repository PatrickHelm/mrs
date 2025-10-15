%neededPI returns the reachable beliefs of being in regime 1 dependend on
%the observed price in period 1

function PI_needed = needed_pi_markov(PI_markov, observed_price,last_price,prices)


t_max=numel(PI_markov);
PI_needed=cell(1,t_max);

past_price=prices==last_price; 
actual_price=prices==observed_price; 

for t=1:t_max
    R1=PI_markov{t}{1};
    R2=PI_markov{t}{2};
    if t==1
        PI_needed{t}=[R1(past_price,actual_price) R2(past_price,actual_price)];
    else
        counter=prices;
        for j=1:t-1
            var=combvec(counter,prices); 
            counter=var;
            var=fliplr(var');
        end
        ind=var(:,1)==last_price&var(:,2)==observed_price ; %First Column is p_t-1 second column is p_t (observed price)
        PI_neededR1=PI_markov{t}{1}(ind,:);
        PI_neededR2=PI_markov{t}{2}(ind,:);
        PI_neededR1=reshape(PI_neededR1,[],1);
        PI_neededR2=reshape(PI_neededR2,[],1);
        PI_needed{t}=[PI_neededR1 PI_neededR2];
    end
end

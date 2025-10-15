% getMins returns the minimal cost from t+1 to calculate the cost of the
% current period t

function mins =get_mins(PI_needed,prices,t,act_pi, act_price, all_mins, inv)

if t==1 %In the first period all minimum values from period t+1=2 are under examination
    mins=all_mins;
else
    block=numel(prices);   % number of possible prices 
    index=find(PI_needed{t}==act_pi); %returns index of the the value 'actPI' in the vector 'PI_needed'; actPi is belief update 
    ende=index*block;
    start=ende-(block-1);
    future_PI=PI_needed{t+1}(start:ende); %'futurePI' is a vector containing all belief update values arising in t+1 from the actual belief update value 'actPi'
    index_2= prices==act_price; %index_2 is the index of the actual price 'actPrice' in the vector of all prices 'prices'
    PI=future_PI(index_2);%PI is the belief update value in t+1 observing the actual Price 'actPrice'
    all_comb_future=combvec(prices,inv,PI_needed{t+1}')'; %Combination of all possible states in t+
    index_3=all_comb_future(:,3)==PI; %index_3 is a logical array which contains a 1 if the computed beleif update 'PI' is equal to the belief update in the vector 'allCombFuture(:,3)'
    mins=all_mins(index_3); %returns the coresponding minimum values
end
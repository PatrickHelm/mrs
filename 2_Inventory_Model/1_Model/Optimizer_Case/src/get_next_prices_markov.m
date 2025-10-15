%get_next_prices computes the probabilities of the next price in t+1 given
%the actual price and the previous price 

function res = get_next_prices_markov(next_price, p, observed_price,last_price,prices, t)

if t>=1
    var=prices;
    for i=1:t
        var=combvec(var, prices);
    end
    var=fliplr(var'); %Creation of all combinations of P1, P2, ..., PN
    pInd1=var(:,1)==last_price & var(:,2)==observed_price;
   res=next_price(pInd1,:);
else %t=1
    for i=1:numel(prices)
        if last_price==prices(i)
            res=next_price(i,:);
            break
        end
    end
end


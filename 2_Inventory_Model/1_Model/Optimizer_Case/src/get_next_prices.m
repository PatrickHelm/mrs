%get_next_prices computes the probabilities of the next price in t+1 given
%the actual price p

function res = get_next_prices(next_price, p, observed_price, prices, t)

if t>=2
    var=prices;
    for i=1:t-1
        var=combvec(var, prices);
    end
    var=fliplr(var'); %Creation of all combinations of P1, P2, ..., PN
    pInd1=var(:,1)==observed_price;
    res=next_price(pInd1,:);
else %t=1
    for i=1:numel(prices)
        if p==prices(i)
            res=next_price(i,:);
            break
        end
    end
end


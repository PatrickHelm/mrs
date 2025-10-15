% backward_costFunction: Calculates total cost

function cost = backward_cost_function(y,cp,ch,p,d_max,I_max,I,next_price,mins,demand_dist,dist)

d_clean=d_max+1; %0 can also be a observed demand
cost=0;
period_cost=single_period_cost(y,cp,ch,p,d_max,I,demand_dist,dist); % Costs, that are calculated if t would be the last period

% Probability that Inventory gets 0
probIgets0=sum(demand_dist(I+y+1:d_clean));

% Probability that inventory becomes > 0 -> [d=m ... d=2 d=1 d=0]
probIgreater0=0;
i=1;
if I+y <=d_max
    for d=I+y-1:-1:0
        probIgreater0(i)=demand_dist(d+1); % ProbIGreater0 = [d=m ... d=2 d=1 d=0] % +1, because indexing of Matlab
        i=i+1;
    end
else
    probIgreater0=demand_dist; % Probability of highest demand at the beginning
end
prob_price_mat=probIgreater0'*next_price; %Matrix that contains multiplications of probabilities, that demand becomes >0 and the probability of the next Price

% Start with highest demand = d_max
block=numel(next_price); % = Number of different prices
ende=block*(I+y+1); % Last value of all minimum costs from t+1, that are neede for the calculation
i=1;

if I+y==0 % i.e. Inventory in t+1 is always 0 -> probIGets0=1;
    for i=1:ende
        cost=cost+probIgets0*next_price(i)*mins(i);
    end
elseif I+y>0 && I+y <= I_max && I+y-d_max<=0; %I+y can not be >I_max
    for i=1:block
       cost=cost+probIgets0*next_price(i)*mins(i);
    end
    j=i+1;
    for P=1:numel(probIgreater0) %Number of possible scenarios, where the inventory is >0 in t+1
        for np=1:numel(next_price)
            cost=cost+prob_price_mat(P,np)*mins(j);
            j=j+1;
        end
    end  
else % I+y>I_max
    j=(I+y-d_max)*block+1;
    for d=d_clean:-1:1 % d is needed for indexing, so d_clean
        if I+y-(d-1) <=I_max % real value of d -> d_clean-1
            for np=1:numel(next_price)
                cost=cost+prob_price_mat(d,np)*mins(j);
                j=j+1;
            end
        else
            cost=cost+period_cost;
            return;
        end
    end
end

cost=cost+period_cost;
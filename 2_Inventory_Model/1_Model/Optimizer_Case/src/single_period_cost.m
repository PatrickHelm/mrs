%single_period_cost calculates the cost for the last period given an order

function costs = single_period_cost(y,cp,ch,p,d_max,I, demand_dist, dist)

costs=p*y; 

if strcmp(dist, 'uniform') 
    for d=0:d_max
        costs=costs + ch/(d_max+1)*positive_or_zero(I+y-d) + cp/(d_max+1)*positive_or_zero(d-I-y);
    end
else 
    for d=0:d_max
        dProb=demand_prob(demand_dist, d); 
        costs=costs + ch*dProb*positive_or_zero(I+y-d) + cp*dProb*positive_or_zero(d-I-y);
    end
end


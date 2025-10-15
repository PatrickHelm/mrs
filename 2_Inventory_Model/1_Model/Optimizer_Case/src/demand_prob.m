%normProb returns the probability of value val for xy distributed
%values in x

function res = demand_prob(x, val)

res=x(val+1); % +1, because demand is starting at 0